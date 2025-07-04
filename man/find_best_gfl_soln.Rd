% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_best_gfl_soln.R
\name{find_best_gfl_soln}
\alias{find_best_gfl_soln}
\title{Determine the GFL problem which has the best solution}
\usage{
find_best_gfl_soln(mod_sum, penfit = "AICc")
}
\arguments{
\item{mod_sum}{an object outputted from \code{tabulate_solution_path()}.}

\item{penfit}{a string indicating which penalized fit measure will be used to
select the models: one of "AICc" or "BICc".}
}
\value{
A numeric value(s) indicating the index (indicies) of the GFL
problem(s) whose solution contains the model with the best penalized fit.
}
\description{
\code{find_best_gfl_soln()} returns the index/indicies of the GFL problem for
which the solution contains the the model with the best penalized fit overall.
}
