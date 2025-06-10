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
#' @param trt_codes a data frame containing the original treatment
#'   codes (`code`) and names (`trt`). This is required if the original data set
#'   contains the treatment names only.
#'
#' @return a data frame that contains data in `gfl_dat`, along with appended
#' grouped parameter codes.
#' @export
#'
recode_gfl_params <- function(param_est, gfl_dat, par_groups, trt_codes) {

  par_config <- get_groupings(param_est, par_groups, pool_format = FALSE)

  # Set up new treatment columns
  gfl_dat$treat1_gfl <- NA
  gfl_dat$treat2_gfl <- NA

  # Extract treatment groupings
  par_config_d <- par_config[grep('d', par_config)]
  for(i in 1:length(par_config_d)) {
    par_config_d[[i]] <- as.numeric(gsub("d_", "", par_config_d[[i]]))
  }

  #
  if(is.character(gfl_dat$treat1)) {
    gfl_dat$treat1 <- factor(gfl_dat$treat1, levels = trt_codes[, "trt"])
    gfl_dat$treat2 <- factor(gfl_dat$treat2, levels = trt_codes[, "trt"])
  }

  # Re-code treatments
  for(i in 1:dim(gfl_dat)[1]) {
    for(j in 1:length(par_config_d)) {
      if(as.numeric(gfl_dat$treat1[i]) %in% par_config_d[[j]]) {
        gfl_dat$treat1_gfl[i] <- j
      }
      if(as.numeric(gfl_dat$treat2[i]) %in% par_config_d[[j]]) {
        gfl_dat$treat2_gfl[i] <- j
      }
    }
  }

  if(length(grep("b", names(mod$par_groups))) > 0) {

    # Set up new class columns
    for(beta in 1:(length(par_groups) - 1)) {
      gfl_dat[i, paste0("treat1_beta", beta, "_gfl")] <- NA
      gfl_dat[i, paste0("treat3_beta", beta, "_gfl")] <- NA
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

  }

  gfl_dat

}
