#' Recode parameters based on configuration suggested by selected GFL model
#'
#' @description
#' `recode_gfl_params()` returns a matrix containing codes for the pooled
#' parameter groups as determined by a selected GFL-NMR model.
#'
#' @param param_est a vector containing the parameter estimates of the desired
#'   pooled model.
#' @param gfl_dat data frame from `prep_gfl_data`.
#' @param par_groups a list of parameter indices (excluding reference) in the
#'   full model.
#'
#' @return a data frame that contains data in `gfl_dat`, along with appended
#' grouped parameter codes.
#' @export
#'
recode_gfl_params <- function(param_est, gfl_dat, par_groups) {

  par_config <- get_groupings(param_est, par_groups, pool_format = FALSE)

  # Extract treatment groupings
  par_config_d <- par_config[grep('d', par_config)]
  for(i in 1:length(par_config_d)) {
    par_config_d[[i]] <- as.numeric(gsub("d_", "", par_config_d[[i]]))
  }

  # Re-code treatments
  for(i in 1:dim(gfl_dat)[1]) {
    for(j in 1:length(par_config_d)) {
      if(gfl_dat$treat1[i] %in% par_config_d[[j]]) {
        gfl_dat$treat1_gfl[i] <- j
        }
      if(gfl_dat$treat2[i] %in% par_config_d[[j]]) {
        gfl_dat$treat2_gfl[i] <- j
        }
    }
  }

  # Extract covariate-interaction groupings for covariates
  par_config_beta <- list()
  for(beta in 1:(length(par_groups) - 1)) {
    par_config_beta[[beta]] <- par_config[grep(paste0('b', beta), par_config)]
    for(j in 1:length(par_config_beta[[beta]])) {
      par_config_beta[[beta]][[j]] <- c(as.numeric(gsub(paste0('b', beta, '_'),
                                                        "",
                                                        par_config_beta[[beta]][[j]])))
    }
  }

  # Create coding for each covariate's groupings
  for(i in 1:dim(gfl_dat)[1]){
    for(beta in 1:(length(par_groups) - 1)) {
      for(j in 1:length(par_config_beta[[beta]])) {
        if(gfl_dat$treat1[i] %in% par_config_beta[[beta]][[j]]) {
          gfl_dat[i, paste0("treat1_beta", beta, "_gfl")] <- j
          }
        if(gfl_dat$treat2[i] %in% par_config_beta[[beta]][[j]]) {
          gfl_dat[i, paste0("treat2_beta", beta, "_gfl")] <- j
          }
      }
    }
  }

  gfl_dat

}
