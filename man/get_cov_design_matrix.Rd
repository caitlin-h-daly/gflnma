% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cov_design_matrix.R
\name{get_cov_design_matrix}
\alias{get_cov_design_matrix}
\title{Determine the design matrix corresponding to the RTE-covariate interaction
terms}
\usage{
get_cov_design_matrix(
  x_cov,
  class1,
  class2,
  ref,
  center = TRUE,
  rescale = TRUE
)
}
\arguments{
\item{x_cov}{a matrix containing the observed study-level covariate values
for P covariates, where the number of rows is equal to the length of the
inputted \code{y} in \code{solve_gflnma()} and the number of columns is equal to P,
the number of covariates.}

\item{class1}{a numeric (or character) vector with the same length as the
inputted \code{y} in \code{solve_gflnma()}, which contains the codes (or names) of
the comparator's class for each relative treatment effect in \code{y}
(\code{treatment2} vs. \code{treatment1}).}

\item{class2}{a numeric (or character) vector with the same length as the
inputted \code{y} in \code{solve_gflnma()}, which contains the codes (or names) of
the intervention's class for each relative treatment effect in \code{y}
(\code{treatment2} vs. \code{treatment1}).}

\item{ref}{a single character string indicating the code (or name) of the
network reference class.}

\item{center}{either a logical value or a numeric vector  of equal length to
the number of covariates indicating whether covariates should be centered
by their means (TRUE, the default), centered by specific values supplied in
the numeric vector, or not at all (FALSE).}

\item{rescale}{either a logical value or a numeric vector  of equal length to
the number of covariates indicating whether covariates should be scaled
by their standard deviation (TRUE, the default), scaled by specific values
supplied in the numeric vector, or not at all (FALSE).}
}
\value{
A matrix where the number of rows equals the that of \code{x_cov}, and the
number of columns equals (C-1)*P in a network of C unique classes (where C =
T, the number of unique treatments, if interaction terms at treatment-level).
}
\description{
\code{get_cov_design_matrix()} returns a design matrix corresponding to the
relative treatment effect (RTE)-covariate interaction terms at the desired
treatment- (or class-) level. The design matrix is scaled by the standardized
observed covariates.
}
