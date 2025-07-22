
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
#'
#' @return A data frame of occupancy times by state and age.
#' @export
#' @examples
#' Create p_list the tidy way
#'p_df <- transitions_example |>
#'  tidyr::pivot_wider(names_from = from_to, values_from = p) |>
#'  dplyr::arrange(age)

#'age <- p_df$age
#'p_list <- as.list(p_df[ , -1])

# Use existing helpers
#'trans_info <- split_transitions(p_list)
#'from_states <- trans_info$from
#'to_states <- trans_info$to

# Fill and infer
#'p_list <- fill_missing_transitions(p_list, from_states, to_states) |>
#' infer_self_transitions(p_list, from_states, to_states)

# calculate occupancies
#'out <- calc_occupancies(p_list, age, origin_state = "P")
#'head(out)


calc_occupancies <- function(p_list,
                             age,
                             init = NULL,
                             origin_state = NULL,
                             delim = "->",
                             age_interval = 1) {
  n                <- length(age)

  trans_info       <- split_transitions(p_list, delim)
  from_states      <- trans_info$from
  to_states        <- trans_info$to
  transient_states <- from_states
  absorbing_states <- setdiff(to_states, from_states)
  all_states       <- union(from_states, to_states)
  n_trans          <- length(transient_states)
  # Initial distribution
  if (!is.null(init)) {
    stopifnot(all(names(init) %in% transient_states))
    stopifnot(all(init >= 0))
    stopifnot(sum(init) > 0)
    init                     <- init / sum(init)
    init_probs               <- rep(0, n_trans)
    names(init_probs)        <- transient_states
    init_probs[names(init)]  <- init
  } else if (!is.null(origin_state)) {
    stopifnot(origin_state %in% transient_states)
    init_probs               <- rep(0, n_trans)
    names(init_probs)        <- transient_states
    init_probs[origin_state] <- 1
  } else {
    stop("Must supply either `origin_state` or `init`.")
  }

  # Call Rcpp calculator
  state_occup <- calc_occupancy_cpp(
    p_list = p_list,
    transient_states = transient_states,
    all_states = all_states,
    init_probs = init_probs,
    n = n,
    delim = delim
  )

  # scale up/down depending on time step size
  state_occup <- state_occup * age_interval

  # Add age column and return as data.frame
  out <- as.data.frame(state_occup)
  out$age <- age
  out
}


#' Fill in missing transitions with zero probabilities
#'
#' @description Ensures that all possible transitions from a set of `from_states` to a set of `to_states` are present in the `p_list`, filling any missing transitions with zero probability vectors. This omits self-transitions.
#'
#' @param p_list A named list of numeric vectors. Each element corresponds to a transition
#'   between states, e.g., `"P->W"`, `"W->D"`.
#' @param from_states Character vector of origin (transient) states.
#' @param to_states Character vector of destination states (can include absorbing).
#' @param delim Character string separating origin and destination states in the names
#'   of transitions. Default is `"->"`.
#'
#' @return A completed `p_list` list with any missing transitions added as zero vectors.
#' @export

fill_missing_transitions <- function(p_list,
                                     from_states,
                                     to_states,
                                     delim = "->") {
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
#' @param delim Character delimiter used in transition names. Default is `"->"`.
#'
#' @return A modified `p_list` with self-transition vectors (e.g., `"P->P"`) added or updated.
#' @export

infer_self_transitions <- function(p_list,
                                   from_states,
                                   to_states,
                                   delim = "->") {
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
#' @param delim Character delimiter between `from` and `to` in transition names. Default is `"->"`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{`from`}{Character vector of unique origin states}
#'     \item{`to`}{Character vector of unique destination states}
#'   }
#' @export

split_transitions <- function(p_list, delim = "->") {
  trans_names <- names(p_list)
  split_mat   <- do.call(rbind, strsplit(trans_names, delim, fixed = TRUE))
  from_states <- unique(split_mat[, 1])
  to_states   <- unique(split_mat[, 2])
  list(from = from_states, to = to_states)
}
