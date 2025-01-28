#' Produce a table summarizing critical points along one or more solution paths
#'
#' @description
#' `tabulate_solution_path()`...
#'
#' @param x a list of one or more objects of class "genlasso".
#' @param lambda_scale a number indicating the desired rescaling power for the
#'   presentation of lambda (i.e., lambdaRescaled = lambda ^ lambda_scale).
#'   Default is 1.
#' @param lambda_critical a character string, one of "fused", "all", or "best",
#'   indicating whether vertical dashed lines should be drawn at lambdas for
#'   which treatments are fused (or unfused), all hitting and leaving points on
#'   the dual path boundary (see Tibshirani and Taylor 2016), or the best
#'   fitting model based on AICc (i.e., the one with the smallest AICc). Default
#'   is "fuse".
#' @param full_rss a number indicating the RSS of the full model (in the case of
#'   a connected network).
#' @param penfit a string indicating which penalized fit measure will be used to
#'   select the models: one of "AIC", "AICc" (default), "BIC", "BICc".
#' @param penfit_threshold a number specifying the degrees of freedom threshold
#'   for which the penalized fit measure should not be computed. See details.
#' @param digits a vector of numbers indicating desired number of digits for
#'   lambda, rss, and penalized fit measures, respectively.
#' @param par_groups a list of lists of grouped indices corresponding to sets of
#'   of parameters (e.g., d, b1, b2, ...).
#'
#' @return A list or data frame summarizing the degrees of freedom, lambda, RSS,
#' AIC, AICc, BIC, BICc, and parameter groupings.
#' @export
#'
tabulate_solution_path <- function(x,
                                   lambda_scale = 1,
                                   lambda_critical = "fused",
                                   full_rss = NA,
                                   penfit = "AICc",
                                   penfit_threshold = NULL,
                                   digits = c(3, 2, 2),
                                   par_groups = NULL
                                   ) {

  mod_sum_list <- lapply(x,
                    tabulate_single_solution_path,
                    lambda_scale,
                    lambda_critical,
                    full_rss,
                    penfit,
                    penfit_threshold,
                    digits,
                    par_groups)

    mod_sum_df <- do.call(rbind, mod_sum_list)

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
