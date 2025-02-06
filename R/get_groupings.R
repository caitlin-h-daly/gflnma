#' Determine the pooled parameter groupings for a given lambda
#'
#' @description
#' `get_groupings()` groups pooled parameters if their estimates are (nearly)
#' equal, as determined with `all.equal()`.
#'
#' @param par_est a vector of (T-1) + P*(T-1) parameter estimates (relative
#'   effects vs. reference treatment + treatment-covariate interactions) in a
#'   network of T treatments, where P covariates have been adjusted for.
#' @param par_groups a list of lists of grouped indicies corresponding to sets
#'   of parameters (excluding reference) to be pooled within but not between.
#' @param pool_format a logical indicating whether pooled parameters should be
#'   presented as a character string (TRUE) or list (FALSE). Default is TRUE.
#'
#' @return A list of the parameter poolings and the labels of parameters
#' belonging to each group.
#' @export
#'
get_groupings <- function(par_est, par_groups = NULL, pool_format = TRUE) {

  if(is.null(par_groups)) {
    # If list of parameter groups is not provided, create a list of length 1
    par_groups <- list()
    par_groups$Parameters <- 1:length(par_est)
  }

  # Create a nested list of the parameter poolings
  par_pool_list <- list()
  pool_index <- 1

  # To label the parameters properly, keep track of the index of the previously
  # plotted parameters
  par_ind <- 0

  for(i in 1:length(par_groups)) {

    # matrix of parameter indices, estimates, and match indicators
    par_unpooled <- cbind(c(1 + par_ind, par_groups[[i]] + 1),
                          # add 0 for reference parameter
                          c(0, par_est[par_groups[[i]]]),
                          # placeholder for matches indicators
                          FALSE)
    colnames(par_unpooled) <- c("par_index", "par_est", "match")

    while(dim(par_unpooled)[1] > 0) {

      # Calculate differences between all parameter estimates and reference parameter
      diff <- par_unpooled[, "par_est"] - par_unpooled[1, "par_est"]

      # Check to see if the parameter estimates are nearly equal to reference parameter
      par_unpooled[, "match"] <- ifelse(unlist(lapply(diff, all.equal, 0)) == TRUE,
                                        TRUE,
                                        FALSE)

      # Group the nearly equal parameter estimates
      par_pool_list[[pool_index]] <- par_unpooled[which(par_unpooled[, "match"] == TRUE),
                                                  "par_index"]

      # Remaining parameters to group
      par_unpooled <- par_unpooled[-which(par_unpooled[, "par_index"] %in%
                                            par_pool_list[[pool_index]]),
                                   ,
                                   drop = FALSE]

      # Improve parameter labels
      if(names(par_groups)[i] != "Parameters"){
        par_pool_list[[pool_index]] <- paste0(names(par_groups)[i],
                                              "_",
                                              par_pool_list[[pool_index]] - par_ind)
      }

      # Update pooled group index
      pool_index <- pool_index + 1

    }

    # To label the parameters properly, keep track of the index of the
    # previously plotted parameters
    par_ind <- par_ind + length(unlist(par_groups[[i]]))

  }

  # For aesthetics, present pooled parameters in a single character string
  par_pool_str <- list()
  for(i in 1:length(par_pool_list)) {
    par_pool_str[[i]] <- paste0("{",
                                paste0(par_pool_list[[i]], collapse = ","),
                                "}")
  }
  par_pool_str <- paste0(par_pool_str, collapse = ",")

  # Return desired output
  if(pool_format) {
    return(par_pool_str)
  } else {
    par_pool_list
  }

}
