% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_monte_carlo.R
\name{compute_v_approxs}
\alias{compute_v_approxs}
\title{Variance approximation of estimators of total}
\usage{
compute_v_approxs(population, v_approx_fn_list, x_names)
}
\arguments{
\item{population}{A data frame corresponding to the population.}

\item{v_approx_fn_list}{A named list of sampling approximations to apply
to the population.}

\item{x_names}{A vector of the names of the auxiliary variables used to
balance the sample.}
}
\value{
A named vector with the variance approximations. The names
correspond to the names in the list.
}
\description{
Multiple variance approximations can be applied on the population.
}
\examples{
population <- data.frame(y = 1:4, x1 = 5:8, pi_i = rep(0.2, 4))
compute_v_approxs(
  population,
  list(v_dt = compute_v_approx_deville_tille),
  "x1"
)
}
