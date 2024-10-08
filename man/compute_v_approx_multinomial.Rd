% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{compute_v_approx_multinomial}
\alias{compute_v_approx_multinomial}
\title{Multinomial approximation of variance of a balanced sample}
\usage{
compute_v_approx_multinomial(population, x_names, tol = 1e-06)
}
\arguments{
\item{population}{A data frame corresponding to the population.
It has:
\itemize{
\item \code{y}: the variable of interest.
\item \code{pi_i}: the inclusion probability of the individual.
}}

\item{x_names}{A vector of the names of the auxiliary variables used to
balance the sample.}

\item{tol}{Threshold to consider an eigenvalue as null.}
}
\value{
A variance approximation of the total estimator of \code{y}.
}
\description{
Multinomial approximation of variance of a balanced sample
}
