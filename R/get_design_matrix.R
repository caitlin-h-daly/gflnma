#' Determine the design matrix representing the direct evidence in the network
#'
#' @description
#' `get_design_matrix()` returns a design matrix X corresponding to the
#' contrast-based info.
#'
#' @param treatment1 a numeric (or character) vector with the same length as the
#'   inputted `y` in `solve_gflnma()`, which contains the codes (or names) of
#'   the comparator for each relative treatment effect in `y` (`treatment2` vs.
#'   `treatment1`).
#' @param treatment2 a numeric (or character) vector with the same length as the
#'   inputted `y` in `solve_gflnma()`, which contains the codes (or names) of
#'   the intervention for each relative treatment effect in `y` (`treatment2`
#'   vs. `treatment1`).
#' @param ref a single character string indicating the code (or name) of the
#'   network reference treatment.
#'
#' @return A matrix where the number of rows equals the length of `treatment1`
#' (or `treatment2`), and the number of columns equals (T-1) in a network of T
#' unique treatments.
#' @export
#'
get_design_matrix <- function(treatment1, treatment2, ref) {

  # Recode the treatments with sequential numbers to determine their position
  # in the design matrix
  new_treat_codes <- recode_treatments(treatment1, treatment2, ref)

  # Create the design matrix
  X <- matrix(0,
              nrow = nrow(new_treat_codes),
              ncol = length(unique(c(new_treat_codes))))
  for (i in 1:nrow(X)) {
    X[i, new_treat_codes[i, ]] <- c(-1, 1)
  }
  colnames(X) <- 1:length(unique(c(new_treat_codes)))

  # Remove the unnecessary column corresponding to the reference treatment in X
  X <- X[, -1]

  as.matrix(X)

}
