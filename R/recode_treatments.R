#' Recode the treatments in sequential order
#'
#' @description
#' `recode_treatments()` returns a matrix containing new codes for the
#' inputted `treatment1` and `treatment2` which correspond to some observed
#' (relative) treatment effects. The new codes ensure the reference treatment
#' (`ref`) is coded as treatment 1, and the remaining treatments are coded from
#' 2 to T in numerical (or alphabetical) order of their original codes (or
#' names), where T is the number of treatments in the network.
#'
#' @param treatment1 a numeric (or character) vector with the same length as
#'   `treatment2`, which contains the codes (or names) of the comparator for
#'   each observed relative treatment effect (`treatment2` vs. `treatment1`).
#' @param treatment2 a numeric (or character) vector with the same length as
#'   `treatment1`, which contains the codes (or names) of the intervention for
#'   each observed relative treatment effect (`treatment2` vs. `treatment1`).
#' @param ref a single character string indicating the code (or name) of the
#'   network reference treatment.
#'
#' @return A matrix of `length(treatment1)` rows and two columns.
#'
recode_treatments <- function(treatment1, treatment2, ref) {

  # Obtain the unique treatment codes (or names) to determine number of
  # treatments and hence appropriate (sequential) treatment codes
  treat_levels <- sort(unique(c(treatment1, treatment2)))
  treat_levels <- c(ref, treat_levels[-which(treat_levels == ref)])  # To code reference treatment as 1

  # Reorder treatment1 and treatment2 and assign them their unique code
  treat1 <- factor(treatment1,
                   levels = treat_levels,
                   labels = seq(1, length(treat_levels), 1))
  treat2 <- factor(treatment2,
                   levels = treat_levels,
                   labels = seq(1, length(treat_levels), 1))

  cbind(treat1, treat2)

}
