#' Determine the GFL problem which has the best solution
#'
#' @description
#' `find_best_gfl()` returns the index/indicies of the GFL problem for which the
#' solution contains the the model with the best penalized fit overall.
#'
#' @param mod_sum_list a list of all model summaries obtained from
#'   `tabulate_solution_path()`.
#' @param mod_sum_df a data frame of all model summaries obtained from
#'   `tabulate_solution_path()`.
#' @param penfit a string indicating which penalized fit measure will be used to
#'   select the models: one of "AICc" or "BICc".
#'
#' @return A numeric value(s) indicating the index (indicies) of the GFL
#' problem(s) whose solution contains the model with the best penalized fit.
#' @export
#'
find_best_gfl <- function(mod_sum_list, mod_sum_df, penfit = "AICc") {

  # Only mod_sum_df has recalculated delta`penfit` statistics across all GFL
  # problems (mod_sum_list calculates delta`penfit` within GFL problems), so we
  # first need to obtain best GFL index from this object
  best_gfl_df_index <- which(mod_sum_df[, paste0("delta", penfit)] == 0)

  # Now find the model(s) in mod_sum_list has the same `penfit` as the one that
  # has the best delta`penfit` in mod_sum_df
  best_gfl_index <- lapply(mod_sum_list,
                             function(x) { which(x[, penfit] ==
                                                   mod_sum_df[, penfit][best_gfl_df_index]) })

  # Extract index/indicies of non-empty sublists
  best_gfl_index <- which(vapply(best_gfl_index,
                                 function(x) { length(x) > 0 }, NA) == TRUE)

  best_gfl_index

}
