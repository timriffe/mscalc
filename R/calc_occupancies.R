
#' Calculate multistate occupancy times
#'
#' @description Calculates state occupancy over age using a list of transition probabilities.
#'
#' @param p_list A named list of age-specific transition probabilities. Recall data.frames and tibbles can function as lists too.
#' @param age A numeric vector of age values.
#' @param init Optional named vector of initial state composition.
#' @param origin_state Optional single origin state (alternative to `init`).
#' @param delim Character delimiter separating from/to in transition names.
#' @param age_interval integer value representing time step (default 1)
#' @param age_col Name of the age column (default = `"age"`) when input is a data.frame.
#' @param trans_col Name of the transition column (default = `"from_to"`) when input is a tidy `data.frame`.
#' @param p_col Name of the transition probability column (default = `"p"`) when input is a tidy `data.frame`.
#' @return A data frame of occupancy times by state and age.
#' @export
#' @importFrom stats setNames
#' @examples
#'
#'prep <- prepare_p_list(transitions_example)
#'occup <- calc_occupancies(prep$p_list, prep$age, origin_state = "P")
#'head(occup)
#' # or without prep (but with prep overhead happening internally)
#' occup <- calc_occupancies(transitions_example, origin_state = "P")
#' # it's ca 2x faster to have the data in list format in advance, FYI.
calc_occupancies <- function(p_list,
                             age = NULL,
                             init = NULL,
                             origin_state = NULL,
                             delim = NULL,
                             age_interval = 1,
                             age_col = "age",
                             trans_col = "from_to",
                             p_col = "p") {


  # Handle inputs flexibly
  is_df <- inherits(p_list, "data.frame")
  is_numeric_list <- is.list(p_list) && !is_df && all(sapply(p_list, is.numeric))

  if (!is_numeric_list) {
    prep <- prepare_p_list(p_list,
                           delim = delim,
                           age_col = age_col,
                           trans_col = trans_col,
                           p_col = p_col)
    age <- prep$age
    p_list <- prep$p_list
  } else if (is.null(age)) {
    stop("If supplying p_list as a list, you must also supply `age`.")
  }


  n <- length(age)

  # Parse transitions
  if (is.null(delim)) {
    delim <- detect_delim(names(p_list))
  }

  trans_info <- split_transitions(p_list, delim)
  from_states <- trans_info$from
  to_states <- trans_info$to
  transient_states <- from_states
  absorbing_states <- setdiff(to_states, from_states)
  all_states <- union(from_states, to_states)

  # Validate / build initial state distribution
  if (!is.null(init)) {
    stopifnot(all(names(init) %in% transient_states))
    stopifnot(abs(sum(init) - 1) < 1e-8)
    init_probs <- setNames(rep(0, length(transient_states)), transient_states)
    init_probs[names(init)] <- init
  } else if (!is.null(origin_state)) {
    stopifnot(origin_state %in% transient_states)
    init_probs <- setNames(rep(0, length(transient_states)), transient_states)
    init_probs[origin_state] <- 1
  } else {
    stop("Must supply either `origin_state` or `init`.")
  }

  # NA check
  all_probs <- unlist(p_list)
  if (any(is.na(all_probs))) {
    stop("Transition probabilities must not contain NA values.")
  }

  # Length check
  lengths_vec <- vapply(p_list, length, integer(1))
  if (length(unique(lengths_vec)) != 1 || unique(lengths_vec) != n) {
    stop("All transition vectors in `p_list` must have the same length as `age`.")
  }

  # Probability bounds check
  if (any(all_probs < 0 | all_probs > 1)) {
    stop("Transition probabilities must be between 0 and 1 (inclusive).")
  }

  # Initialize occupancy matrix
  state_occup <- matrix(0, nrow = n, ncol = length(transient_states))
  colnames(state_occup) <- transient_states
  state_occup[1, ] <- init_probs

  # Main loop
  for (i in 1:(n - 1)) {
    current <- state_occup[i, ]
    next_occup <- setNames(numeric(length(transient_states)), transient_states)

    for (from in transient_states) {
      occup_from <- current[from]
      outflow <- 0

      for (to in all_states) {
        trans_name <- paste(from, to, sep = delim)
        p <- if (trans_name %in% names(p_list)) p_list[[trans_name]][i] else 0
        if (to %in% transient_states) {
          next_occup[to] <- next_occup[to] + occup_from * p
        }
        outflow <- outflow + p
      }

      # Residual self-transition
      stay_prob <- max(0, 1 - outflow)
      next_occup[from] <- next_occup[from] + occup_from * stay_prob
    }

    state_occup[i + 1, ] <- next_occup
  }

  # Adjust based on age interval
  state_occup <- state_occup * age_interval

  out <- as.data.frame(state_occup)
  out$age <- age
  return(out)
}
#' Fill in missing transitions with zero probabilities
#'
#' @description Ensures that all possible transitions from a set of `from_states` to a set of `to_states` are present in the `p_list`, filling any missing transitions with zero probability vectors. This omits self-transitions.
#'
#' @param p_list A named list of numeric vectors. Each element corresponds to a transition
#'   between states, e.g., `"P->W"`, `"W->D"`.
#' @param from_states Character vector of origin (transient) states.
#' @param to_states Character vector of destination states (can include absorbing).
#' @param delim Character string separating origin and destination states in the names of transitions. If not specified we try to detect it from names of `p_list`
#'
#' @return A completed `p_list` list with any missing transitions added as zero vectors.
#' @export

fill_missing_transitions <- function(p_list,
                                     from_states,
                                     to_states,
                                     delim = NULL) {
  if (is.null(delim)){
    delim <- detect_delim(names(p_list))
  }
  n <- length(p_list[[1]])
  all_pairs <- expand.grid(from = from_states,
                           to = to_states,
                           stringsAsFactors = FALSE)
  for (i in seq_len(nrow(all_pairs))) {
    name <- paste(all_pairs$from[i], all_pairs$to[i], sep = delim)
    if (!name %in% names(p_list)) {
      p_list[[name]] <- rep(0, n)
    }
  }
  p_list
}


#' Infer self-transitions as residual probabilities
#'
#' @description For each transient state, calculates the remaining probability after all defined outflows and assigns it to the self-transition (e.g., `"P->P"`).
#'
#' @param p_list A named list of numeric vectors representing transition probabilities.
#' @param from_states Character vector of origin states (transient).
#' @param to_states Character vector of possible destination states.
#' @param delim Character delimiter used in transition names. If not specified we try to detect it from names of `p_list`
#'
#' @return A modified `p_list` with self-transition vectors (e.g., `"P->P"`) added or updated.
#' @export

infer_self_transitions <- function(p_list,
                                   from_states,
                                   to_states,
                                   delim = NULL) {
  if (is.null(delim)){
    delim <- detect_delim(names(p_list))
  }
  n <- length(p_list[[1]])
  all_states <- union(from_states, to_states)
  for (from in from_states) {
    outflows <- rep(0, n)
    for (to in all_states) {
      if (to != from) {
        name <- paste(from, to, sep = delim)
        if (name %in% names(p_list)) {
          outflows <- outflows + p_list[[name]]
        }
      }
    }
    self_name <- paste(from, from, sep = delim)
    p_list[[self_name]] <- pmax(0, 1 - outflows)
  }
  p_list
}
#' Extract origin and destination states from transition names
#'
#' @description Parses the names of transitions in `p_list` to identify the unique sets of origin (`from`) and destination (`to`) states.
#'
#' @param p_list A named list of transitions where names follow the `"from->to"` convention.
#' @param delim Character delimiter between `from` and `to` in transition names. If not specified we try to detect it from names of `p_list`
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`from`}{Character vector of unique origin states}
#'     \item{`to`}{Character vector of unique destination states}
#'   }
#' @export

split_transitions <- function(p_list, delim = NULL) {
  if (is.null(delim)){
    delim <- detect_delim(names(p_list))
  }
  trans_names <- names(p_list)
  split_mat   <- do.call(rbind, strsplit(trans_names, delim, fixed = TRUE))
  from_states <- unique(split_mat[, 1])
  to_states   <- unique(split_mat[, 2])
  list(from = from_states, to = to_states)
}

#' Prepare transition probability list for occupancy calculations
#'
#' Prepares a list of transition probability vectors (`p_list`) suitable for input to
#' multistate occupancy functions. Accepts tidy long-format data, wide-format data,
#' or a named list. Automatically fills in missing transitions and infers self-transitions.
#'
#' @param transitions Either:
#'   - A tidy data frame with columns `age`, `from_to`, and `p`,
#'   - A wide data frame with `age` and one column per transition (e.g. "PW", "WD"),
#'   - Or a named list of transition vectors, with names of the form "from_to".
#' @param age_col The name of the age column (default is `"age"`).
#' @param trans_col The name of the column containing information of origin-destination pairs for transitions. Default `"from_to"`.
#' @param p_col The name of the column containing transition probabilities. Default `"p"`.
#' @param delim Character string used to split transition names into origin and destination states. Default is `"->"`, but use `""` if transitions are compact like `"PW"` or `"UD"`.
#'
#' @return A list with two elements:
#'   - `p_list`: A named list of transition probabilities with all required transitions included.
#'   - `age`: A numeric vector of ages, extracted from input data if applicable.
#'
#' @examples
#' data(transitions_example)
#'
#' # From tidy input
#' prep <- prepare_p_list(transitions_example, delim = "")
#'
#' # From wide data frame
#' wide <- tidyr::pivot_wider(transitions_example,
#'  names_from = from_to,
#'  values_from = p)
#' prep2 <- prepare_p_list(wide, delim = "")
#'
#' # From list input
#' age <- wide$age
#' p_list <- as.list(wide[, -1])
#' prep3 <- prepare_p_list(p_list, delim = "")
#'
#' @export

prepare_p_list <- function(transitions,
                           age_col = "age",
                           trans_col = "from_to",
                           p_col = "p",
                           delim = NULL) {

  # Infer delimiter if not provided
  if (is.null(delim)) {
    name_sample <- if (is.data.frame(transitions)) {
      if (trans_col %in% names(transitions)) {
        unique(transitions[[trans_col]])
      } else {
        names(transitions)[names(transitions) != age_col]
      }
    } else if (is.list(transitions)) {
      names(transitions)
    } else {
      stop("Cannot detect delimiter: unrecognized input format.")
    }
    delim <- detect_delim(name_sample)
  }

  # Case 1: It's already a named list of numeric vectors
  if (is.list(transitions) && !is.data.frame(transitions)) {
    p_list <- transitions
    age <- seq_along(transitions[[1]])
  }
  # Case 2: Tidy data frame (long format)
  else if (is.data.frame(transitions) &&
           all(c(age_col, trans_col, p_col) %in% names(transitions))) {

    wide <- tidyr::pivot_wider(transitions,
                               id_cols = all_of(age_col),
                               names_from = all_of(trans_col),
                               values_from = all_of(p_col)) |>
      dplyr::arrange(.data[[age_col]])

    age <- wide[[age_col]]
    p_list <- as.list(wide[ , !(names(wide) %in% age_col)])
  }
  # Case 3: Wide data frame with age + transition columns + optional metadata
  else if (is.data.frame(transitions) && age_col %in% names(transitions)) {
    wide <- dplyr::arrange(transitions, .data[[age_col]])
    age <- wide[[age_col]]

    # Filter only transition columns
    maybe_trans <- setdiff(names(wide), age_col)
    if (length(maybe_trans) == 0) stop("No transition columns detected in wide format.")

    if (is.null(delim)) {
      delim <- detect_delim(maybe_trans)
    }

    # Only include column names that match the state transition pattern
    if (delim == "") {
      # No delimiter: only allow names of length 2 (e.g., "PW")
      trans_cols <- maybe_trans[nchar(maybe_trans) == 2]
    } else {
      transition_pattern <- paste0("^[^", delim, "]+", delim, "[^", delim, "]+$")
      trans_cols <- grep(transition_pattern, maybe_trans, value = TRUE)
    }

    if (length(trans_cols) == 0) {
      stop("No valid transition columns detected in wide format.")
    }


    p_list <- as.list(wide[ , trans_cols, drop = FALSE])
  } else {
    stop("Unsupported input format for `transitions`. Must be tidy df, wide df, or list.")
  }

  # Transition name check
  trans_info <- split_transitions(p_list, delim = delim)
  from_states <- trans_info$from
  to_states <- trans_info$to

  # Degenerate case: only self-transitions (e.g., PP)
  if (all(from_states == to_states)) {
    if (delim == "") delim <- "->"

    # Convert to standard form with absorbing state
    new_p_list <- list()
    for (name in names(p_list)) {
      from <- substr(name, 1, nchar(name) / 2)
      p_self <- p_list[[name]]
      new_p_list[[paste(from, from, sep = delim)]] <- p_self
      new_p_list[[paste(from, "abs", sep = delim)]] <- pmax(0, 1 - p_self)
    }
    p_list <- new_p_list

    trans_info <- split_transitions(p_list, delim)
    from_states <- trans_info$from
    to_states   <- trans_info$to
  }

  # Complete transition structure
  p_list <- fill_missing_transitions(p_list, from_states, to_states, delim)
  p_list <- infer_self_transitions(p_list, from_states, to_states, delim)

  return(list(p_list = p_list, age = age))
}

#' Detect delimiter used in transition names
#'
#' Attempts to infer the delimiter used in the names of transition probability vectors.
#'
#' @param names A character vector of transition names (e.g. `"P->W"`, `"PW"`, `"P_W"`).
#'
#' @return A character string representing the detected delimiter, or `""` if none found.
#'
#' @examples
#' detect_delim(c("P->W", "W->D"))     # returns "->"
#' detect_delim(c("PW", "WR", "RD"))   # returns ""
#' detect_delim(c("P_W", "W_D"))       # returns "_"
#'
#' @export
detect_delim <- function(names) {
  if (!is.character(names)) stop("Input to detect_delim() must be a character vector.")

  common_delims <- c("->", "_", "-", ":", "::", "\\|", "~", "\\.")

  detected <- vapply(common_delims, function(d) {
    all(grepl(d, names))
  }, logical(1))

  if (any(detected)) {
    return(common_delims[which(detected)[1]])
  } else {
    return("")  # Assume compact form like "PW" with no delimiter
  }
}

