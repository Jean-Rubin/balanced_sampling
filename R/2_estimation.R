#' Estimate total from a sample
#'
#' The sample can have non-integer selection of units, for example as a result
#' of the flight phase.
#'
#' @param sample A data frame corresponding to a sample of the population.
#'   It contains:
#'   - `y`: the variable of interest.
#'   - `pi_i`: the inclusion probability.
#'   - `s_i`: the (eventually non-integer) indicator of inclusion.
#'
#' @return A scalar corresponding to the estimator of the total.
#' @examples
#' samp <- data.frame(y = 1:5, pi_i = rep(0.5, 5), s_i = c(1, 0, 1, 0.5, 0, 0))
#' estimate_total(samp)
estimate_total <- function(sample) {
  sample |>
    summarize(y_hat = sum(s_i * y / pi_i)) |>
    pull(y_hat)
}

estimate_variance_balanced <- function(sample_pop, x_names, tol = 1e-6) {
  x <- as.matrix(sample_pop[x_names])
  y <- as.matrix(sample_pop$y)
  pi_i <- sample_pop$pi_i
  x_diag_p <- t(x / pi_i)
  beta_hat <- pseudo_inv(x_diag_p %*% x, tol) %*% x_diag_p %*% y
  e <- (y - x %*% beta_hat) / pi_i

  t(e) %*% e
}
