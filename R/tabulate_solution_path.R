#' Produce a table summarizing critical points along one or more solution paths
#'
#' @description
#' `tabulate_solution_path()`...
#'
#' @param x a list of one or more objects of class "genlasso".
#' @param lambda_scale a number indicating the desired rescaling power for the
#'   presentation of lambda (i.e., lambdaRescaled = lambda ^ lambda_scale).
#'   Default is 1.
#' @param lambda_scale_digits a number indicating desired number of digits for
#'   rescaled lambda. Default is 3.
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
#'   penalized fit measures. Default is 4.
#' @param par_groups a list of lists of grouped indices corresponding to sets of
#'   of parameters (e.g., d, b1, b2, ...).
#'
#' @return A list or data frame summarizing the degrees of freedom, lambda, RSS,
#' AIC, AICc, BIC, BICc, and parameter groupings.
#' @export
#'
tabulate_solution_path <- function(x,
                                   lambda_scale = 1,
                                   lambda_scale_digits = 3,
                                   lambda_critical = "fused",
                                   full_rss = NA,
                                   penfit = "AICc",
                                   penfit_threshold = NULL,
                                   penfit_digits = 4,
                                   par_groups = NULL
                                   ) {

  mod_sum_list <- lapply(x,
                    tabulate_single_solution_path,
                    lambda_scale,
                    lambda_scale_digits,
                    lambda_critical,
                    full_rss,
                    penfit,
                    penfit_threshold,
                    penfit_digits,
                    par_groups)

    mod_sum_df <- do.call(rbind, mod_sum_list)

    # If the network is connected (via x$bls the least square solution), add
    # summary statistics for full model
    if(!is.null(x[[1]]$bls) & is.na(full_rss)) {
      warning("RSS of full model not specified. Please specify to include full
            model summary statistics in table.")
    } else if(!is.null(x[[1]]$bls) & !is.na(full_rss)) {
      full_df <- dim(x[[1]]$beta)[1]
      full_rss <- as.numeric(full_rss)
      full_AIC <- (2 * full_df) + full_rss
      full_AICc <- full_AIC + (2 * full_df * (full_df + 1)) /
        (length(x[[1]]$y) - full_df - 1)
      full_BIC <- (log(length(x[[1]]$y)) * full_df) + full_rss
      full_BICc <- full_BIC + (2 * full_df * (full_df + 1)) /
        (length(x[[1]]$y) - full_df - 1)
      full_param_groups <- paste0("{",
                                  rep(names(par_groups),
                                      each = length(par_groups$d) + 1),
                                  "_",
                                  1:(length(par_groups$d) + 1),
                                  "}",
                                  collapse = ",")
      mod_sum_df[nrow(mod_sum_df) + 1,] <- c(NA,  # gamma
                                             NA,  # eps
                                             full_df,  # degrees of freedom
                                             0,  # lambda
                                             full_rss,  # RSS
                                             full_AIC,  # AIC
                                             NA,  # deltaAIC
                                             full_AICc,  # AICc
                                             NA,  # deltaAICc
                                             full_BIC,  # BIC
                                             NA,  # deltaBIC
                                             full_BICc,  # BICc
                                             NA,  # deltaBICc
                                             NA  # Pooled groups
                                             )
      mod_sum_df[nrow(mod_sum_df), "PooledGroups"] <- full_param_groups
    }

    # Re-calculate delta penalty statistics over all gammas and eps
    mod_sum_df[, "deltaAIC"] <- mod_sum_df[, "AIC"] -
      mod_sum_df[, "AIC"][which(mod_sum_df[, "AIC"] == min(mod_sum_df[, "AIC"]))]
    mod_sum_df[, "deltaAICc"] <- mod_sum_df[, "AICc"] -
      mod_sum_df[, "AICc"][which(mod_sum_df[, "AICc"] == min(mod_sum_df[, "AICc"]))]
    mod_sum_df[, "deltaBIC"] <- mod_sum_df[, "BIC"] -
      mod_sum_df[, "BIC"][which(mod_sum_df[, "BIC"] == min(mod_sum_df[, "BIC"]))]
    mod_sum_df[, "deltaBICc"] <- mod_sum_df[, "BICc"] -
      mod_sum_df[, "BICc"][which(mod_sum_df[, "BICc"] == min(mod_sum_df[, "BICc"]))]

  list("mod_sum_list" = mod_sum_list, "mod_sum_df" = mod_sum_df)

}
