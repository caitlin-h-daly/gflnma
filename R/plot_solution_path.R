#' Plot the primal solution path
#'
#' @description
#' `plot_solution_path()` returns a plot of the primal solution path for objects
#' of class `genlasso`. User must provide the `genlasso` object obtained from
#' `solve_gflnma()`.
#'
#' @details
#' To draw vertical dashed lines at critical lambda values for which the
#' parameters pool or unpool, the user must provide a `genlasso_summary` object
#' obtained from `tabulate_solution_path()` via the `gfl_sum` argument. The user
#' may draw additional lines at specific lambda values on an ad-hoc basis via
#' `abline()`.
#'
#' The user also has the option to specify a critical region of `penfit` values
#' through `delta_penfit_crit`. See Kong et al. (2024) for recommended rules
#' of thumb. Additionally, the user has the option to specify which sets of
#' parameters should be plotted separately by specifying a listed of their
#' grouped indices through `par_groups`.
#'
#' @param gfl_soln an object of class `genlasso`.
#' @param mod_index an optional numeric value indicating the index of the
#'   GFL-NMA problem to plot, if multiple problems were solved by
#'   `solve_gflnma()` (e.g.,if `length(gamma)` > 1).
#' @param par_groups an optional list of lists of grouped indices corresponding
#'   to sets of parameters (excluding reference) to be plotted separately
#'   (e.g., d, b1, b2, ...).
#' @param gfl_sum an optional object of class `genlasso_summary` produced by the
#'   `tabulate_solution_path()` function, which is used to determine the
#'   critical lambda values for which parameters pool or unpool.
#' @param penfit an optional string indicating which penalized fit measure will
#'   be used to select the models: one of "AIC", "AICc", "BIC", "BICc".
#' @param delta_penfit_lim an optional number indicating the range of
#'   delta_`penfit` values corresponding to models the user is willing to
#'   consider. Default is NULL.
#' @param lambda_scale a number indicating the desired rescaling power for the
#'   presentation of lambda (i.e., lambda_rescaled = lambda ^ `lambda_scale`).
#'   Default is 1.
#' @param xlimits numeric vector of length 2, giving the x coordinates range.
#' @param ylimits numeric vector of length 2, giving the y coordinates range.
#' @param cex_lab a numeric value indicating the expansion factor for the axis
#'   titles.
#' @param cex_axis a numeric value indicating the expansion factor for the axis
#'   labels.
#' @param cex_text a numeric value indicating the expansion factor for the
#'   superimposed parameter labels.
#'
#' @return One or more plots of the estimated parameters against lambda (the
#' penalty factor for d).
#' @export
#'
plot_solution_path <- function(gfl_soln,
                               mod_index = NULL,
                               par_groups = NULL,
                               gfl_sum = NULL,
                               penfit = NULL,
                               delta_penfit_lim = NULL,
                               lambda_scale = 1,
                               xlimits = NULL,
                               ylimits = NULL,
                               cex_lab = 2,
                               cex_axis = 2,
                               cex_text = 0.8) {

  if (length(gfl_soln) > 1) {
    stop("Multiple GFL-NMA problems detected. Please specify which problem you
         would like to plot the solution for via `mod_index`.")
  }

  if(!is.null(mod_index)) {
    gfl_soln <- gfl_soln[[mod_index]]
    if(!is.null(gfl_sum)) {
      gfl_sum <- gfl_sum[[mod_index]]
    }
  } else {
    gfl_soln <- gfl_soln[[1]]
    if(!is.null(gfl_sum)) {
      gfl_sum <- gfl_sum[[1]]
    }
  }

  if (!inherits(gfl_soln, "genlasso"))
    stop("\'gfl_soln\' must be an object of class `genlasso`.")

  if (!is.null(gfl_sum) & !inherits(gfl_sum, "genlasso_summary"))
    stop("\'gfl_sum\' must be an object of class `genlasso_summary` obtained
         from `genlasso_summary_table()`.")

  # Create string for x-axis label if lambda is rescaled
  xlabel <- if(lambda_scale != 1) {
    bquote(lambda  ^ .(lambda_scale))
    } else {
      bquote(lambda)
    }

  # Set up plotting range for x-axis
  if(is.null(xlimits)) {
    xlimits <- c(0, gfl_soln$lambda[1] * 1.1)
  } else {
    xlimits <- xlimits
  }

  if(is.null(par_groups)) {
    # If list of parameter groups is not provided, create a list of length 1
    par_groups <- list()
    par_groups$Parameters <- 1:nrow(gfl_soln$beta)
  }

  # To label the parameters properly, keep track of the index of the previously
  # plotted parameters
  par_ind <- 0

  # Now plot the solution paths
  for(i in 1:length(par_groups)) {

    # Line colours
    cols <- rainbow(n = nrow(gfl_soln$beta[par_groups[[i]], ]),
                    start = 1/6)  # exclude red

    # If GLS solutions exists, include this in the plot
    if(exists("bls", gfl_soln)) {
      x_pts <- c(gfl_soln$lambda, 0)
      y_pts <- rbind(t(gfl_soln$beta[par_groups[[i]], ]),
                     t(gfl_soln$bls[par_groups[[i]]]))
    } else {
      x_pts <- gfl_soln$lambda
      y_pts <- t(gfl_soln$beta[par_groups[[i]], ])
    }

    # Set up plotting range for y-axis
    if(is.null(ylimits)) {
      ylimits_temp <- range(y_pts) + c(-0.1, 0.1)
    } else {
      ylimits_temp <- ylimits
    }

    par(mar = c(5.1, 5.1, 4.1, 2.1))

    # Solution path
    matplot(x = x_pts,
            y = y_pts,
            type = "l",
            col = cols,
            lty = 1,
            xlab = xlabel,
            ylab = paste0("Coordinates of ", names(par_groups)[i]),
            xlim = xlimits,
            ylim = ylimits_temp,
            cex.lab = cex_lab,
            cex.axis = cex_axis)

    # Add horizontal line for reference parameter (d_1 = b1_1 = b2_1 = ... = 0)
    abline(h = 0, col = "red")

    if(!is.null(penfit) & is.null(gfl_sum)){
      warning("`gfl_sum` must be specified to plot vertical line at the lambda
              with smallest penalized fit.")
    }

    if(!is.null(gfl_sum)) {

      # Add vertical dashed lines at every lambda where parameters pool or unpool
      abline(v = gfl_sum$lambda ^ lambda_scale, lty = 2)

      if(!is.null(penfit)) {
        # Add vertical red line at lambda with smallest penfit
        best_lambda_index <- which(gfl_sum[, paste0("delta", penfit)] == 0)
        abline(v = gfl_sum$lambda[best_lambda_index] ^ lambda_scale,
               col = "red",
               lty = 2)
      }

    }

    # Add parameter labels to lines
    if(names(par_groups)[i] == "Parameters") {
      text(x = xlimits[1],
           y = y_pts[nrow(y_pts),],
           labels = par_groups[[i]] + 1,
           cex = cex_text)
    } else {
      text(x = xlimits[1],
           y = y_pts[nrow(y_pts),],
           labels = paste0(names(par_groups)[i], "_", (par_groups[[i]] + 1) - par_ind),
           cex = cex_text)
    }

    legend("right",
           # List treatments in order of parameter estimates at smallest lambda
           legend = paste0(names(par_groups)[i],
                           "_",
                           order(y_pts[nrow(y_pts),], decreasing = TRUE) + 1),
           # Assign colours by parameter estimates at smallest lambda
           col = cols[order(y_pts[nrow(y_pts),], decreasing = TRUE)],
           lty = 1)

    # Optional: Shade the region containing the lambda values that correspond to
    # models with delta_penfit within the critical region
    if(!is.null(gfl_sum) & !is.null(delta_penfit_lim)) {
      # Index of models that fall within delta_penfit_lim range
      penfit_lim_index <- which(gfl_sum[, paste0("delta", penfit)] > 0 &
                                  gfl_sum[, paste0("delta", penfit)] < abs(delta_penfit_lim))
      rect(xleft = min(gfl_sum$lambda[penfit_lim_index]) ^ lambda_scale,
           xright = max(gfl_sum$lambda[penfit_lim_index]) ^ lambda_scale,
           ybottom = min(y_pts[nrow(y_pts),]) - diff(range(y_pts[nrow(y_pts),])),
           ytop = max(y_pts[nrow(y_pts),]) + diff(range(y_pts[nrow(y_pts),])),
           border = NA,
           col = adjustcolor("blue", alpha = 0.2))

    }

    # To label the parameters properly, keep track of the index of the
    # previously plotted parameters
    par_ind <- par_ind + length(par_groups[[i]])

  }

}
