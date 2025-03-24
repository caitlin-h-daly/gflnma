#' Determine the design matrix representing the direct evidence in the network
#'
#' @description
#' `get_cov_design_matrix()` returns a design matrix corresponding to the
#' treatment- (or class-) covariate interaction terms at the desired treatment-
#' or class-level. The design matrix is scaled by the standardized
#' observed covariates.
#'
#' @param x_cov a matrix containing the observed study-level covariate values
#'    for P covariates, where the number of rows is equal to the length of the
#'    inputted `y` in `solve_gflnma()` and the number of columns is equal to P,
#'    the number of covariates.
#' @param class1 a numeric (or character) vector with the same length as the
#'   inputted `y` in `solve_gflnma()`, which contains the codes (or names) of
#'   the comparator's class for each relative treatment effect in `y`
#'   (`treatment2` vs. `treatment1`).
#' @param class2 a numeric (or character) vector with the same length as the
#'   inputted `y` in `solve_gflnma()`, which contains the codes (or names) of
#'   the intervention's class for each relative treatment effect in `y`
#'   (`treatment2` vs. `treatment1`).
#' @param ref a single character string indicating the code (or name) of the
#'   network reference class.
#' @param center either a logical value or a numeric vector  of equal length to
#'   the number of covariates indicating whether covariates should be centered
#'   by their means (TRUE, the default), centered by specific values supplied in
#'   the numeric vector, or not at all (FALSE).
#' @param rescale either a logical value or a numeric vector  of equal length to
#'   the number of covariates indicating whether covariates should be scaled
#'   by their standard deviation (TRUE, the default), scaled by specific values
#'   supplied in the numeric vector, or not at all (FALSE).
#'
#' @return A matrix where the number of rows equals the that of `x_cov`, and the
#' number of columns equals (C-1)*P in a network of C unique classes (where C =
#' T, the number of unique treatments, if interaction terms at treatment-level).
#' @export
#'
get_cov_design_matrix <- function(x_cov, class1, class2, ref,
                                  center = TRUE, rescale = TRUE) {

  # Design matrix for contrasts at desired class-level (which could = treatment)
  X_class <- get_design_matrix(class1, class2, ref)
  # X_class <- matrix(replicate(ncol(x_cov), X_class),
  #                   ncol = ncol(x_cov) * ncol(X_class))

  # Center the covariates
  if(center == TRUE) {
    x_cov <- sweep(x_cov, MARGIN = 2, STATS = colMeans(x_cov), FUN = "-")
  } else if(is.numeric(center) & length(center) == ncol(x_cov)) {
    x_cov <- sweep(x_cov, MARGIN = 2, STATS = center, FUN = "-")
  } else if(is.numeric(center) & length(center) != ncol(x_cov)){
    error("length(`center`) does not equal number of columns of `x_cov`")
  }

  # Scale the covariates
  if(rescale == TRUE) {
    x_cov <- sweep(x_cov, MARGIN = 2, STATS = apply(x_cov, 2, sd), FUN = "/")
  } else if(is.numeric(rescale) & length(rescale) == ncol(x_cov)) {
    x_cov <- sweep(x_cov, MARGIN = 2, STATS = rescale, FUN = "/")
  } else if(is.numeric(rescale) & length(rescale) != ncol(x_cov)){
    error("length(`rescale`) does not equal number of columns of `x_cov`")
  }

  # Design matrix that incorporates the covariates to adjust for
  X_cov <- matrix(nrow = nrow(x_cov), ncol = ncol(X_class) * ncol(x_cov))
  for(i in 1:dim(x_cov)[1]) {
    param_index <- 1
    for(p in 1:dim(x_cov)[2]) {
      X_cov[i, param_index:(param_index + ncol(X_class) - 1)] <- x_cov[i,p] * X_class[i,]
      param_index <- param_index + ncol(X_class)
    }
  }

  colnames(X_cov) <- paste0("beta",
                            rep(1:dim(x_cov)[2], each = ncol(X_class)),
                            "_",
                            1:ncol(X_class) + 1)

  X_cov

}
