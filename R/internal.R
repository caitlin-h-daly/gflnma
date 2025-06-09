#' Produce a table summarizing critical points along a single solution path
#'
#' @description
#' `tabulate_single_solution_path` calls on `genlasso::summary.genlasso` to
#' tabulate the degrees of freedom (df), the tuning parameter for the relative
#' treatment effects (lambda), and residual sum of squares (RSS), the
#' corresponding penalized model fit statistics (AIC, AICc, BIC, BICc), and the
#' parameter groupings at each critical lambda where parameters pool or
#' unpool.
#'
#' @details
#' In star networks with one study per comparison, the degrees of freedom (df)
#' equals the number of studies (M) when no treatments are pooled. The AIC and
#' BIC measures corrected for small sample size are minimized when df = M, and
#' so the full model will always be selected in these networks. Additionally,
#' the correction term for small sample size in AIC and BIC will be undefined
#' when df = M - 1. Thus, for these networks, we recommend not considering these
#' models. By default, the table will exclude these models in the summary
#' output. Alternatively, the user can specify `penfit_threshold` such that the
#' table will exclude models in the summary output for which
#' (df+1)/M >= `penfit_threshold`. If the user wants all models to be presented,
#' set `penfit_threshold` to be any number > (df+1)/M.
#'
#' @param x an object of class "`genlasso`".
#' @param lambda_scale a number indicating the desired rescaling power for the
#'   presentation of lambda (i.e., lambdaRescaled = lambda ^ lambda_scale).
#'   Default is 1.
#' @param lambda_scale_digits a number indicating desired number of digits for
#'   rescaled lambda.
#' @param lambda_critical a character string, one of "fused", "all", or "best",
#'   indicating whether vertical dashed lines should be drawn at lambdas for
#'   which treatments are fused (or unfused), all hitting and leaving points on
#'   the dual path boundary (see Tibshirani and Taylor 2016), or the best
#'   fitting model based on AICc (i.e., the one with the smallest AICc). Default
#'   is "fuse".
#' @param full_rss a number indicating the RSS of the full model (if estimable).
#' @param penfit a string indicating which penalized fit measure will be used to
#'   select the models: one of "AIC", "AICc" (default), "BIC", "BICc".
#' @param penfit_threshold a number specifying the degrees of freedom threshold
#'   for which the penalized fit measure should not be computed. See details.
#' @param penfit_digits a number indicating desired number of digits for the
#'   penalized fit measures.
#' @param par_groups a list of lists of grouped indices corresponding to sets of
#'   of parameters (excluding reference) to be plotted separately.
#'
#' @return A data frame summarizing the df, lambda, RSS, AIC, AICc and treatment
#'   groupings.
#'
tabulate_single_solution_path <- function(x,
                                          lambda_scale,
                                          lambda_scale_digits,
                                          lambda_critical,
                                          full_rss,
                                          penfit,
                                          penfit_threshold,
                                          penfit_digits,
                                          par_groups
                                          ) {

  if(list(NULL) %in% x) {
    return(NULL)
  } else {

    if (!inherits(x, "genlasso"))
      stop("\'x\' must include an object of class 'genlasso', or an object which
         inherits this class (i.e., 'fusedlasso').")

    # Obtain summary statistics which include required RSS
    fit_summary <- genlasso::summary.genlasso(x)

    # Start building the table through a list
    ll <- list()

    # Denote gamma and eps values
    if(exists(x = "gamma", where = x) & exists(x = "eps", where = x)) {
      ll[["gamma"]] <- x[["gamma"]]
      ll[["eps"]] <- x[["eps"]]
    } else {
      ll[["gamma"]] <- NA
      ll[["eps"]] <- NA
    }

    # Extract statistics required to calculate AICc
    ll[["df"]] <- fit_summary[, "df"]
    ll[["lambda"]] <- fit_summary[, "lambda"]
    ll[["rss"]] <- fit_summary[, "rss"]

    # Remove any models that do not meet the threshold for consideration
    if(is.null(penfit_threshold)) {
      # If the user has not specified a threshold, check to see if df >= M - 1
      # (which may be the case in star networks)
      if(any((ll[["df"]] + 1) / length(x$y) >= 1)) {
        warning("Note models with df = M or df = M + 1 (M = number of studies) are
              not presented. See Details for more information.")
        # Note index of models that should be removed
        ind_penfit_threshold <- which((ll[["df"]] + 1) / length(x$y) >= 1)
        # Remove these models from the fit statistics
        ll <- lapply(ll, function(x) x[-ind_penfit_threshold])
        # Remove the coefficient estimates from these models as well
        #if(!is.null(x$bls) & !is.na(full_rss)){
        # if the statistics for the full model have been added, it should not be
        # indexed for removal from x$beta
        #x[["beta"]] <- x[["beta"]][, -ind_penfit_threshold[-length(ind_penfit_threshold)]]
        #} else {
        x[["beta"]] <- x[["beta"]][, -ind_penfit_threshold]
        #}
      }
    } else {
      # Check to see if any of the models fall outside the user's set consider-
      # ation range
      if(any((ll[["df"]] + 1) / length(x$y) >= penfit_threshold)) {
        warning(paste("Note models with (df+1)/M >=", penfit_threshold,
                      "are not presented. See Details for more information."))
        ind_penfit_threshold <- which((ll[["df"]] + 1) / length(x$y) >= penfit_threshold)
        # Remove these models from the fit statistics
        ll <- lapply(ll, function(x) x[-ind_penfit_threshold])
        # Remove the coefficient estimates from these models as well
        #if(!is.null(x$bls) & !is.na(full_rss)) {
        # if the statistics for the full model have been added, it should not be
        # indexed for removal from x$beta
        #  x[["beta"]] <- x[["beta"]][, -ind_penfit_threshold[-length(ind_penfit_threshold)]]
        # } else {
        x[["beta"]] <- x[["beta"]][, -ind_penfit_threshold]
        #}
      }
    }

    # Calculate AIC and deltaAIC
    ll[["AIC"]] <- (2 * ll[["df"]]) + ll[["rss"]]
    ll[["deltaAIC"]] <- ll[["AIC"]] - ll[["AIC"]][which(ll[["AIC"]] == min(ll[["AIC"]]))]

    # Calculate AICc and deltaAICc
    ll[["AICc"]] <- ll[["AIC"]] + ((2 * ll[["df"]] * (ll[["df"]] + 1)) / (length(x$y) - ll[["df"]] - 1))
    ll[["deltaAICc"]] <- ll[["AICc"]] - ll[["AICc"]][which(ll[["AICc"]] == min(ll[["AICc"]]))]

    # Calculate BIC and deltaBIC
    ll[["BIC"]] <- (log(length(x$y)) * ll[["df"]]) + ll[["rss"]]
    ll[["deltaBIC"]] <- ll[["BIC"]] - ll[["BIC"]][which(ll[["BIC"]] == min(ll[["BIC"]]))]

    # Calculate BICc and deltaBICc
    ll[["BICc"]] <- ll[["BIC"]] + ((2 * ll[["df"]] * (ll[["df"]] + 1)) / (length(x$y) - ll[["df"]] - 1))
    ll[["deltaBICc"]] <- ll[["BICc"]] - ll[["BICc"]][which(ll[["BICc"]] == min(ll[["BICc"]]))]

    # Add lists of pooled parameters
    ll[["PooledGroups"]] <- unname(apply(x$beta, 2, FUN = get_groupings, par_groups))

    # Rescale lambda
    if(lambda_scale != 1) {
      ll[[paste0("lambda_pow",lambda_scale)]] <- ll[["lambda"]] ^ lambda_scale
      if(!is.null(lambda_scale_digits)) {
        ll[[paste0("lambda_pow",lambda_scale)]] <- round(ll[[paste0("lambda_pow", lambda_scale)]],
                                                         digits = lambda_scale_digits)
      }
    }

    # Round numeric values
    if(!is.null(penfit_digits)) {
      ll[["AIC"]] <- round(ll[["AIC"]], digits = penfit_digits)
      ll[["AICc"]] <- round(ll[["AICc"]], digits = penfit_digits)
      ll[["deltaAIC"]] <- round(ll[["deltaAIC"]], digits = penfit_digits)
      ll[["deltaAICc"]] <- round(ll[["deltaAICc"]], digits = penfit_digits)
      ll[["BIC"]] <- round(ll[["BIC"]], digits = penfit_digits)
      ll[["BICc"]] <- round(ll[["BICc"]], digits = penfit_digits)
      ll[["deltaBIC"]] <- round(ll[["deltaBIC"]], digits = penfit_digits)
      ll[["deltaBICc"]] <- round(ll[["deltaBICc"]], digits = penfit_digits)
    }

    # Convert list into data.frame
    tt <- do.call(data.frame, ll)

    # Temporary indicator function to group dfs for which the solution path has
    # the same slope
    ## First ensure tt is ordered by descending lambda
    tt <- tt[order(tt[, "lambda"], decreasing = TRUE), ]
    ## Now create group indicator
    df_group_dum <- rep(NA, dim(tt)[1])
    df_group_dum[1] <- 1
    for(i in 2:dim(tt)[1])
      df_group_dum[i] <- ifelse(tt[i, "df"] == tt[i-1, "df"],
                                df_group_dum[i-1],
                                df_group_dum[i-1] + 1)
    ## Rows that contain smallest lambda for each df group
    df_group_min <- NA
    for(i in 1:max(df_group_dum)) {
      df_group_min[i] <- max(which(df_group_dum == i))
    }

    # Subset rows based on desired output indicated by lambda_critical
    if(lambda_critical == "fused") {
      tt <- tt[df_group_min, ]
    } else if(lambda_critical == "best") {
      tt <- tt[which(tt[, penfit] == min(tt[, penfit])), ]
    } else {
      tt <- tt
    }

    # Sort table by increasing lambda
    tt <- tt[order(tt[, "lambda"], decreasing = FALSE), ]

    # Specify class of tt which will be checked in the genlasso_path_plot()
    # function
    class(tt) <- c("genlasso_summary", "data.frame")

    # remove row names
    row.names(tt) <- NULL

    tt

  }

}


