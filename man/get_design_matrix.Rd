% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_design_matrix.R
\name{get_design_matrix}
\alias{get_design_matrix}
\title{Determine the design matrix representing the direct evidence in the network}
\usage{
get_design_matrix(treatment1, treatment2, ref)
}
\arguments{
\item{treatment1}{a numeric (or character) vector with the same length as the
inputted \code{y} in \code{solve_gflnma()}, which contains the codes (or names) of
the comparator for each relative treatment effect in \code{y} (\code{treatment2} vs.
\code{treatment1}).}

\item{treatment2}{a numeric (or character) vector with the same length as the
inputted \code{y} in \code{solve_gflnma()}, which contains the codes (or names) of
the intervention for each relative treatment effect in \code{y} (\code{treatment2}
vs. \code{treatment1}).}

\item{ref}{a single character string indicating the code (or name) of the
network reference treatment.}
}
\value{
A matrix where the number of rows equals the length of \code{treatment1}
(or \code{treatment2}), and the number of columns equals (T-1) in a network of T
unique treatments.
}
\description{
\code{get_design_matrix()} returns a design matrix X corresponding to the
contrast-based info.
}
