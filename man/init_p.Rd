% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonlinear_optimisation.R
\name{init_p}
\alias{init_p}
\title{Produce an initial probability vector for the max-entropy sampling}
\usage{
init_p(lambda, eps = 1e-11)
}
\arguments{
\item{eps}{Threshold to consider a value as zero.}
}
\value{
A vector of probability with the same dimension as lambda.
}
\description{
Every probability is initialized with a non-zero value.
}
\examples{
init_p(c(1, 2, -3))
}
