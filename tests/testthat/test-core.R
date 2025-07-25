test_that("Transition probabilities sum to 1 after filling self-transitions", {
  p_list <- prepare_p_list(transitions_example)
  trans_info <- split_transitions(p_list$p_list)
  from_states <- trans_info$from
  to_states <- trans_info$to
  n <- length(p_list$p_list[[1]])

  for (from in from_states) {
    total_outflow <- rep(0, n)
    for (to in union(from_states, to_states)) {
      trans_name <- paste0(from, to)
      if (trans_name %in% names(p_list$p_list)) {
        total_outflow <- total_outflow + p_list$p_list[[trans_name]]

      }
    }
    expect_true(all(abs(total_outflow - 1) < 1e-8))
  }
})


test_that("Negative probabilities trigger error", {
  bad_list <- prepare_p_list(transitions_example)
  bad_list$p_list[["PW"]][10] <- -0.01
  expect_error(calc_occupancies(bad_list$p_list,
                                origin_state = "P",
                                age = bad_list$age))
})


test_that("Delimiter detection works for different styles", {
  trans_names <- c("P->D", "P->W", "P->R")
  expect_equal(detect_delim(trans_names), "->")

  trans_names <- c("P_W", "P_R", "P_D")
  expect_equal(detect_delim(trans_names), "_")

  trans_names <- c("P.R", "P.D", "P.W")
  expect_equal(detect_delim(trans_names), "\\.")

  trans_names <- c("PW", "PR", "PD")
  expect_equal(detect_delim(trans_names), "")
})


test_that("Unlinked destination states don't cause failure", {
  test_list <- list(PW = rep(0.5, 10), PD = rep(0.5, 10), ZZ = rep(0, 10))
  expect_silent(split_transitions(test_list))
  expect_silent(fill_missing_transitions(test_list, c("P"), c("W", "D", "Z")))
})


test_that("calc_occupancies requires either init or origin_state", {
  p_list <- prepare_p_list(transitions_example)
  age <- 65:100
  expect_error(calc_occupancies(p_list, age = age))
})

test_that("1-state with PP or PD works like lifetable", {
  age <- 65:70
  # P -> P (perfect survivorship)
  p_list <- list(PP = rep(.9, length(age)))
  p_list <- prepare_p_list(transitions = p_list,age=age)
  res <- calc_occupancies(p_list$p_list, age = age, origin_state = "P", delim="->")
  expect_equal(res$P, c(1,cumprod(rep(.9, length(age)-1))))

  # P -> D (some death)
  p_list <- list(PD = rep(0.1, length(age)))
  p_list <- prepare_p_list(transitions = p_list,age=age)
  res2 <- calc_occupancies(p_list$p_list, age = age, origin_state = "P")
  expect_true(all(diff(res2$P) < 0))  # should decline
})

test_that("Single-state system with no transitions throws error", {
  p_list <- list()  # no transitions
  age <- 65:70
  expect_error(calc_occupancies(p_list, age = age, origin_state = "P"))
})

test_that("NA transition probabilities trigger error", {
  bad_list <- prepare_p_list(transitions_example)
  bad_list$p_list[["PW"]][5] <- NA
  expect_error(calc_occupancies(bad_list$p_list,
                                origin_state = "P",
                                age = bad_list$age),
               regexp = "NA values")
})

test_that("Mismatched transition length triggers error", {
  bad_list <- prepare_p_list(transitions_example)
  bad_list$p_list[["PW"]] <- bad_list$p_list[["PW"]][-1]  # shorten vector
  expect_error(calc_occupancies(bad_list$p_list,
                                origin_state = "P",
                                age = bad_list$age),
               regexp = "length")
})

test_that("Unlinked states do not cause failure", {
  p_list <- list(PW = rep(0.5, 6), WX = rep(0.3, 6))  # 'X' not used
  age <- 65:70
  out <- calc_occupancies(p_list, age = age, origin_state = "P")
  expect_true(is.data.frame(out))
})

test_that("Tidy input is accepted and produces expected results", {
  result <- calc_occupancies(transitions_example, origin_state = "P")
  expect_s3_class(result, "data.frame")
  expect_true("P" %in% names(result))
  expect_true("age" %in% names(result))
})

test_that("Wide data frame input is handled", {
  wide <- tidyr::pivot_wider(transitions_example,
                             names_from = from_to,
                             values_from = p)
  wide <- dplyr::arrange(wide, age)

  result <- calc_occupancies(wide, origin_state = "P")
  expect_s3_class(result, "data.frame")
  expect_true("P" %in% names(result))
})


test_that("Tidy input with renamed columns is handled", {
  df <- transitions_example |>
    dplyr::rename(age_group = age,
                  trans_label = from_to,
                  probability = p)

  out <- calc_occupancies(df,
                          origin_state = "P",
                          age_col = "age_group",
                          trans_col = "trans_label",
                          p_col = "probability")

  expect_true("P" %in% colnames(out))
  expect_equal(nrow(out), length(unique(df$age_group)))
})

test_that("Wide input with renamed age column and extra metadata is handled", {
  wide <- tidyr::pivot_wider(transitions_example,
                             names_from = from_to,
                             values_from = p)
  wide <- dplyr::rename(wide, start_age = age)
  wide$region <- "X"  # extra irrelevant metadata

  out <- calc_occupancies(wide,
                          origin_state = "P",
                          age_col = "start_age")

  expect_true("P" %in% colnames(out))
  expect_equal(nrow(out), length(unique(wide$start_age)))
})

