create_pseudo_population <- function(sample_pop, x_names) {
  n_sample <- nrow(sample_pop)
  sample_pop |>
    tidyr::uncount(.env$n_sample) |>
    mutate(
      across(all_of(x_names), \(x) x / pi_i),
      y = y / (.env$n_sample * pi_i),
      pi_i = 1 / .env$n_sample
    )
}


## Calibration approach -----
# with-replacement sampling with a multinomial reweighting
# calibrate on t_x with x the auxiliary variable
calibrate_up <- function(sample_pop, bal_var) {
  X <- as.matrix(sample_pop[bal_var])
  dim_X <- dim(X)
  if ((dim_X[1] == 0) || (dim_X[2] == 0)) return(sample_pop)
  sample_pop <- sample_pop |>
    mutate(
      st_i = 1,
      G_i = n() / (n() - 1) * rmultinom(1, n() - 1, rep(1 / n(), n())),
      d_i = 1 / pi_i
    )

  d_star <- sample_pop$G_i * sample_pop$d_i
  nb_st <- table(sample_pop$st_i)
  for (h in names(nb_st)) {
    in_h <- sample_pop$st_i == h
    X_h <- X[in_h, , drop = FALSE]
    t_z <- t(X_h) %*% sample_pop$d_i[in_h]
    t_z_star <- t(X_h) %*% d_star[in_h]
    A_hat <- t(X_h) %*% (X_h * d_star[in_h])
    lambda <- solve(A_hat) %*% (t_z - t_z_star)
    sample_pop$G_i[in_h] <- sample_pop$G_i[in_h] * (1 + X[in_h, ] %*% lambda)
  }

  sample_pop |>
    mutate(
      pi_i = 1 / (as.numeric(G_i) * d_i),
      pi_i_aux = pi_i
    )
}
