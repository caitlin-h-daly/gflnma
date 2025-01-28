#' Determine the penalty matrix to implement GFL within parameter groups
#'
#' @description
#' `get_penalty_matrix()` determines a matrix of the penalty terms for the
#' relative treatment effects and the treatment-covariate interaction terms (if
#' applicable), where the latter is scaled by the desired ratio of penalization
#' factors (`gamma`, see details for more information).
#'
#' @details
#' Currently designed for GFL-NMA where the relative treatment effects against a
#' reference treatmetn (d) and their pairwise differences are penalized, as well
#' as GFL-NMR which additionally (but separately) penalizes the treatment-
#' covariate interactions (b), as well as their pairwise differences.
#'
#' Regarding the ratio of penalty factors: if gamma > 1, then the treatment-
#' covariate interactions will be penalized gamma times more than the treatment
#' effects. If gamma = 0, the treatment-covariate interactions are not
#' penalized. Treatment effects and their pairwise differences receive the same
#' penalty, while treatment-covariate interactions and their pairwise
#' differences receive gamma*that penalty.
#'
#' @param X a design matrix corresponding to the observed relative effects
#'   in `y` inputted into `solve_gflnma()`, which may be obtained from
#'   `get_design_matrix()`. This matrix may also be a binded matrix of the
#'   design matrices for the relative effects and the interaction terms, the
#'   latter of which may be obtained from `get_covariate_design_matrix()`.
#' @param x_cov an optional matrix containing the observed study-level covariate
#'   values for P covariates, where the number of rows is equal to the length of
#'   the inputted `y` in `solve_gflnma()` and the number of columns is equal to
#'   P, the number of covariates. This should be specified if X includes
#'   interaction terms
#' @param gamma a non-negative number indicating the ratio of penalization
#'   factors (treatment-specific covariate interactions (b)/relative treatment
#'   effects (d)). Default is 1.
#'
#' @return A matrix consisting of with (T-1) + P*(T-1) columns, where T is the
#' number of treatments in the network, and P is the number of covariates.
#' @export
#'
get_penalty_matrix <- function(X, x_cov = NULL, gamma = 1) {

  # In case X includes design matrix for interaction terms, extract matrix for
  # `d` terms
  if(!is.null(x_cov)) {
    X <- X[, 1:(ncol(X) / ncol(x_cov))]
  }

  # Create oriented incidence matrix
  penalty_matrix <- matrix(0, nrow = choose(ncol(X), 2), ncol = ncol(X))
  j <- 1
  k <- 1
  for(i in 1:choose(ncol(X), 2)) {
    penalty_matrix[i, c(j, j + k)] <- c(-1, 1)
    k <- k + 1
    if((j + k) > ncol(X)) {
      j <- j + 1
      k <- 1
    }
  }
  # Add column names
  colnames(penalty_matrix) <- paste0("d_", 1:ncol(X) + 1)

  # Append identity matrix to oriented incidence matrix
  penalty_matrix <- rbind(penalty_matrix, diag(ncol(X)))

  if(!is.null(x_cov)) {

    if(gamma < 0)
      stop("`gamma` must be non-negative.")

    # Create oriented matrix for interaction terms
    penalty_matrix <- as.matrix(do.call(Matrix::bdiag,
                                        rep(list(penalty_matrix),
                                            ncol(x_cov))) * gamma)

    # Replace column names
    colnames(penalty_matrix) <- paste0("beta",
                                       rep(1:ncol(x_cov), each = ncol(X)),
                                       "_",
                                       1:ncol(X) + 1)

  }

  penalty_matrix

}
