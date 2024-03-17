estimate_total <- function(sample) {
  sample |>
    summarize(y_hat = sum(y / pi_i * s_i)) |>
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






