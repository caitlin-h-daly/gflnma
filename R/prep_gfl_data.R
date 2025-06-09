#' Prepare contrast-level data for GFL-NMA or GFL-NMR
#'
#' @description
#' `prep_gfl_data()` determines the covariance of relative effects (vs. a common
#' comparator) in multi-arm studies and returns a list of a data frame that only
#' includes the data required by `solve_gflnma()`, as well as the variance-
#' covariance matrix required by `solve_gflnma()`, and a string indicating the
#' assumed type of effect for between-study heterogeneity.
#'
#' @details
#' User must specify if they want to fit a fixed effect (FE), random effects
#' (RE), or multiplicative effects (ME) model. This is critical for
#' determining the correct variance-covariance matrix.
#'
#' `prep_gfl_data()` currently calculates covariance between contrasts within
#' studies using the variance of the outcome in the common comparator group via
#' one of the following:
#'   * `se1`, the standard error of the observed estimates in the comparator
#'   group. If all pairwise relative effects are available, this argument may be
#'   used to input the covariance via equation (4.12) of Dias et al. (2018).
#'   * `sd1` and `n1`, the standard deviation and sample size of the observed
#'   outcomes in the comparator group.
#'   * `event1` and `time1`, the number of observed events and the person time
#'   at risk in the comparator groups.
#'   * `event1` and `n1`, the number of observed events and the number of
#'   participants in the comparator groups.
#'
#' `prep_gfl_data()` currently calculates the random effect of between-study
#' heterogeneity based on a random effects model fitted via
#' `netmeta::netmeta()`. We recommend using multiplicative effects GFL-NMR over
#' random effects GFL-NMR.
#'
#' @references Dias S, Ades AE, Welton NJ, Jansen JP, Sutton AJ. Network meta-
#' analysis for decision-making. John Wiley & Sons; 2018.
#'
#' @param dat a data frame containing the observed data in each study. This must
#'   contain all pairwise comparisons within multiarm studies (e.g., B vs. A,
#'   C vs. A, C vs. B) if a random effects model is desired, so that an estimate
#'   of tau can be obtained from `netmeta::netmeta()`.
#' @param studlab the name of the column in `dat` that contains the study
#'   labels.
#' @param TE the name of the column in `dat` that contains the estimates of the
#'   relative treatment effects.
#' @param seTE the name of the column in `dat` that contains the standard errors
#'   of the relative treatment effects.
#' @param treatment1 the name of the column in `dat` that contains the labels of
#'   the comparator for each relative treatment effect in `TE` (`treatment2`
#'   vs. `treatment1`).
#' @param treatment2 the name of the column in `dat` that contains the labels of
#'   the intervention for each relative treatment effect in `TE` (`treatment2`
#'   vs. `treatment1`).
#' @param V the name of the column in `dat` that contains the covariance of the
#'   of the observed relative effects.
#' @param se1 the name of the column in `dat` that contains the standard error
#'   of the observed estimates in the comparator arm.
#' @param sd1 the name of the column in `dat` that contains the standard
#'   deviation of the observed outcomes in the comparator arm.
#' @param n1 the name of the column in `dat` that contains the number of
#'   participants in the comparator arm.
#' @param event1 the name of the column in `dat` that contains the number of
#'   events observed among participants in the comparator arm.
#' @param time1 the name of the column in `dat` that contains the total person
#'   time at risk among participants in the comparator arm.
#' @param mod_type a character string indicating if a fixed effect (FE), random
#'   effects (RE), or multiplicative effects (ME) model is to be fitted.
#'
#' @return A list of
#'   * a data frame containing the required study-level data for GFL-NMA or
#'     GFL-NMR,
#'   * a C x C covariance matrix for the observed relative effects, where C is
#'     the number of rows in the outputted data frame,
#'   * a character string indicating the model assumed for the covariance
#'     in the outputted data frame.
#' @export
#'
prep_gfl_data <- function(dat,
                          TE,
                          studlab,
                          seTE,
                          treatment1,
                          treatment2,
                          V = NULL,
                          se1 = NULL,
                          sd1 = NULL,
                          n1 = NULL,
                          event1 = NULL,
                          time1 = NULL,
                          modtype = "ME") {

  # Find most frequent comparator (and hence baseline treatment) in each study
  study_comp <- table(dat[, studlab], dat[, treatment1])
  study_base <- data.frame(study = row.names(study_comp),
                           base_trt = colnames(study_comp)[apply(study_comp, 1, which.max)])

  # Denote the baseline treatment in each trial
  for(i in 1:nrow(dat)){
    dat$base_trt[i] <- study_base[which(study_base[, "study"] == dat[i, studlab]), "base_trt"]
  }

  # Keep a copy of the data set with all pairwise comparisons for netmeta
  dat_netmeta <- dat
  # Keep the contrasts where the comparator is the baseline treatment
  dat <- dat[which(dat[, treatment1] == dat$base_trt), ]

  # Sort the data by study label + reference treatment of the contrast
  dat <- dat[order(dat[, studlab], dat[, treatment1]), ]
  row.names(dat) <- NULL

  # Calculate the covariances of the observed contrasts within multi-arm studies
  # For each trial, denote the number of arms (n_arms), contrast index
  # (contr_ind), column index of covariances (cov_ind2), and covariance of
  # relative effects (cov), which will be useful when we populate the var-cov
  # matrix
  dat$n_arms <- NA
  dat$contr_ind <- NA
  dat$contr_ind2 <- NA
  dat$cov <- ifelse(!is.null(V), dat[, V], 0)
  for(i in 1:nrow(dat)) {

    dat_sub <- dat[which(dat[, studlab] == dat[i, studlab]), ]

    dat$n_arms[i] <- nrow(dat_sub) + 1

    # To determine indices of covariance entries for multi-arm studies, need to
    # keep track of index of study-specific contrasts
    if(i == 1) {
      dat$contr_ind[i] <- 1
    } else if(dat[i, studlab] == dat[i - 1, studlab]) {
      dat$contr_ind[i] <- dat$contr_ind[i - 1] + 1
    } else {
      dat$contr_ind[i] <- 1
    }

    # Matrix of indices for off-diagonals
    delta <- matrix(rep(seq_len(nrow(dat_sub)), nrow(dat_sub)) -
                      rep(seq_len(nrow(dat_sub)), each = nrow(dat_sub)),
                    nrow = nrow(dat_sub),
                    byrow = TRUE)

    if(nrow(dat_sub) > 1) {
      dat$contr_ind2[i] <- list(i + delta[dat$contr_ind[i],
                                        which(delta[dat$contr_ind[i], ] != 0)])
    }

    if(!is.null(V)) {
      dat$cov[i] <- ifelse(is.na(dat[i, V]), 0, dat[i, V])
    } else{
      # Add covariance of relative effects, which is equal to the variance of the
      # arm-level estimate in the baseline arm
      if(dat$n_arms[i] > 2 & dat$base_trt[i] == dat[i, treatment1]) {
        if(!is.null(se1)) {
          dat$cov[i] <- dat[i, se1]^2
        } else if(!is.null(sd1) & !is.null(n1)) {
          dat$cov[i] <- (dat[i, sd1]^2)/dat[i, n1]
        } else if(!is.null(event1) & !is.null(time1)) {
          dat$cov[i] <- 1/dat[i, event1]
        } else if(!is.null(event1) & !is.null(n1)) {
          dat$cov[i] <- dat[i, event1]*(dat[i, n1] - dat[i, event1]) / (dat[i, n1] ^ 3)
        } else if(dat$n_arms[i] > 2 & dat$base_trt[i] != dat[i, treatment1]) {
          dat$cov[i] <- "Something went wrong"
        }
      }
    }
  }

  if(modtype == "FE" | modtype == "ME") {

    # Fixed or mixed effects var-cov: matrix with all var on diagonal
    var_cov <- diag((dat[, seTE])^2)
    # Add the covariances for multi-arm trials
    for(i in 1:nrow(dat)) {
      if(dat$n_arms[i] > 2) {
        var_cov[i, unlist(dat$contr_ind2[i])] <- dat$cov[i]
      }
    }

    mod_type <- ifelse(modtype == "FE", "fixed", "multiplicative")
  }

  if(modtype == "RE") {

    # To fit a random effects GFL-NMR model, add tau to the covariance.
    # To obtain an estimate of tau, we need to fit a RE model via `netmeta`.

    # Obtain tau from NMA model fitted via netmeta.
    net <- netmeta::netmeta(TE = dat_netmeta[, TE],
                            seTE = dat_netmeta[, seTE],
                            treat1 = dat_netmeta[, treatment1],
                            treat2 = dat_netmeta[, treatment2],
                            studlab = dat_netmeta[, studlab],
                            common = FALSE)

    # Covariances under RE model.
    dat$cov_RE <- ifelse(dat$n_arms > 2, dat$cov + (net$tau2 / 2), 0)

    # Random effects var-cov: matrix with all var on diagonal
    var_cov <- diag( ((dat[, seTE])^2) + net$tau2 )
    # Add the covariances for multi-arm trials
    for(i in 1:nrow(dat)) {
      if(dat$n_arms[i] > 2) {
        var_cov[i, unlist(dat$contr_ind2[i])] <- dat$cov_RE[i]
      }
    }

    mod_type <- "random"

  }

  # Clean up data frame
  dat$base_trt <- NULL
  dat$n_arms <- NULL
  dat$contr_ind <- NULL
  dat$contr_ind2 <- NULL

  list(dat = dat, var_cov = var_cov, mod_type = mod_type)

}
