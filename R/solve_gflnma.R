#' Solve one or more GFL-NMA or GFL-NMR models
#'
#' @description
#' `solve_gflnma()` constructs one or more penalty matrix based on inputted
#' `gamma` and calls on `genlasso::genlasso()` to compute the solution path of
#' the generalized fused lasso problem(s).
#'
#' @param y a numeric vector containing the observed relative treatment effects
#'   vs. a study-specific common treatment for each study.
#' @param var_cov a block-diagonal matrix of each study's variance-covariance
#'   matrix corresponding to observed relative treatment effects in `y`, where
#'   the number of rows and number of columns both equal the length of the `y`.
#' @param x_cov a matrix containing the observed study-level covariate values
#'    for P covariates, where the number of rows is equal to the length of `y`
#'    and the number of columns is equal to P, the number of covariates.
#' @param treatment1 a numeric (or character) vector with the same length as
#'   `y`, which contains the codes (or names) of the comparator for each
#'   observed relative treatment effect (`treatment2` vs. `treatment1`).
#' @param treatment2 a numeric (or character) vector with the same length as
#'   `y`, which contains the codes (or names) of the intervention for each
#'   observed relative treatment effect (`treatment2` vs. `treatment1`).
#' @param ref a single character string indicating the code (or name) of the
#'   network reference treatment.
#' @param class1 an optional numeric (or character) vector with the same length
#'   as `y`, which contains the codes (or names) of `treatment1`'s class for
#'   each observed relative treatment effect (`class2` vs. `class1`). Only input
#'   if the treatment-covariate terms should be shared within classes.
#' @param class2 an optional numeric (or character) vector with the same length
#'   as `y`, which contains the codes (or names) of `treatment2`'s class for
#'   each observed relative treatment effect (`class2` vs. `class1`). Only input
#'   if the treatment-covariate terms should be shared within classes.
#' @param ref_class a single character string indicating the code (or name) of
#'   the network reference treatment class.  Only input if the treatment-
#'   covariate terms should be shared within classes.
#' @param fit_full a logical value indicating if the full model should be fitted
#'   (TRUE) or not (FALSE, default).
#' @param minlam a numeric variable indicating the value of tuning parameter for
#'   the d's at which the solution path should terminate. Default is 0.
#' @param gamma a vector of numeric values > 0 indicating the desired ratio of
#'   the penalty for interaction terms (beta's) to the penalty for the relative
#'   treatment effects (d's). Default is 0 to only penalize the relative
#'   treatment effects. If only interaction terms are to be penalized, input
#'   `Inf` (no quotes).
#' @param eps a vector of numeric values > 0 indicating the desired
#'   multiplier(s) for the ridge penalty, which will be implemented if there are
#'   some parameters for which there is no direct nor indirect evidence in the
#'   data. Default is 0.0001.
#' @param center either a logical value or a numeric vector  of equal length to
#'   the number of covariates indicating whether covariates should be centered
#'   by their means (TRUE, the default), centered by specific values supplied in
#'   the numeric vector, or not at all (FALSE).
#' @param rescale either a logical value or a numeric vector  of equal length to
#'   the number of covariates indicating whether covariates should be scaled
#'   by their standard deviation (TRUE, the default), scaled by specific values
#'
#' @return A nested list of
#'   * `solution`: a nested list of multiple `genlasso::genlasso` objects.
#'   * `par_groups`: a list of lists of parameter indices (excluding reference).
#' @export
#'
solve_gflnma <- function(y,
                         var_cov,
                         x_cov = NULL,
                         treatment1,
                         treatment2,
                         ref,
                         class1 = NULL,
                         class2 = NULL,
                         ref_class = NULL,
                         fit_full = FALSE,
                         minlam = 0,
                         gamma = 0,
                         eps = 0.0001,
                         center = TRUE,
                         rescale = TRUE) {

  # Obtain weight matrix via upper triangular factor in a Cholesky decomposition
  # of the inverse var_cov matrix
  U <- chol(solve(var_cov))

  # Determine the design matrix for treatment contrasts
  X_d <- get_design_matrix(treatment1, treatment2, ref)
  colnames(X_d) <- paste0("d_", 1:ncol(X_d) + 1)

  # If covariate-adjustment required, append the design matrix for the
  # interaction terms to the design matrix for the treatment contrasts
  if(!is.null(x_cov)) {
    # Determine the class-covariate or treatment-covariate interaction terms
    if(!is.null(class1) & !is.null(class2) & !is.null(ref_class)) {
      X_beta <- get_cov_design_matrix(x_cov, class1, class2, ref_class,
                                      center, rescale)
    } else {
      X_beta <- get_cov_design_matrix(x_cov, treatment1, treatment2, ref,
                                      center, rescale)
    }
    X <- cbind(X_d, X_beta)
  } else {
    X <- X_d
  }

  # If full model requested
  if(fit_full) {
    param <- U %*% X
    mod_full <- lm( (U %*% y) ~ -1 + param )
  }

  # Now solve `length(gamma)` GFL-NMA problems for each eps
  mod <- list()
  index <- 1

  for(i in 1:length(gamma)) {

    # First need to get gamma-specific penalty matrix
    penalty_matrix <- get_penalty_matrix(X_d)

    # if(gamma[i] == Inf | gamma[i] == "Inf"){
    #   gamma[i] <- Inf
    #   warning("The relative effects have not been penalized.")
    #   penalty_matrix <- matrix(0,
    #                            nrow = nrow(penalty_matrix),
    #                            ncol = ncol(penalty_matrix))
    # }

    if(exists("X_beta")) {

      if(gamma[i] == Inf | gamma[i] == "Inf"){
        gamma[i] <- Inf
        warning("The relative effects have not been penalized.")
        penalty_matrix_beta <- get_penalty_matrix(X_beta, x_cov, gamma = 1)
        penalty_matrix <- cbind(matrix(0,
                                       nrow = nrow(penalty_matrix_beta),
                                       ncol = ncol(penalty_matrix)),
                                penalty_matrix_beta)
      } else {
        penalty_matrix_beta <- get_penalty_matrix(X_beta, x_cov, gamma = gamma[i])
        penalty_matrix <- as.matrix(Matrix::bdiag(penalty_matrix,
                                                  penalty_matrix_beta))
      }

      colnames(penalty_matrix) <- c(paste0("d_", 1:ncol(X_d) + 1),
                                    paste0("beta",
                                           rep(1:dim(x_cov)[2],
                                               each = ncol(X_beta) / dim(x_cov)[2]),
                                           "_",
                                           1:(ncol(X_beta) / dim(x_cov)[2]) + 1))
    }


    for(j in 1:length(eps)) {

      mod[[index]] <- tryCatch(genlasso::genlasso(y = U %*% y,
                                         X = U %*% X,
                                         D = penalty_matrix,
                                         minlam = minlam,
                                         eps = eps[j],
                                         svd = TRUE),
                               error = function(e) {
                                 print(paste0("An error happened when computing the solution path: gamma = ", gamma[i], ", eps = ", eps[j]))
                                 return(NA)
                               }
      )

      if(all(is.na(mod[[index]]))){
        mod[[index]] <- NULL
      } else {
        mod[[index]]["gamma"] <- gamma[i]
        mod[[index]]["eps"] <- eps[j]
        index <- index + 1
      }

    }
  }

  # Note indices of d's and beta's in parameter list so we may consider their summaries separately
  par_groups <- list()
  par_groups[["d"]] <- 1:ncol(X_d)
  if(!is.null(x_cov)) {
    for(i in 1:dim(x_cov)[2]) {
      par_groups[[paste0("b", i)]] <- (i * ncol(X_d) + 1):((i + 1) * ncol(X_d))
    }
  }

  if(fit_full){
    list("solution" = mod, "par_groups" = par_groups, "full_RSS" = deviance(mod_full))
  } else {
    list("solution" = mod, "par_groups" = par_groups)
  }

}


