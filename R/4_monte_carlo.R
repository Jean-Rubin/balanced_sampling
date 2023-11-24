mcmc_population <- function(
  population,
  sample_fn,
  n_iter,
  .progress = TRUE
) {
  purrr::map_dbl(
    seq_len(n_iter),
    \(x) population |>
      sample_fn() |>
      estimate_total(),
    .progress = .progress
  )
}

compute_v_trues <- function(population, sample_fn_list, n_iter_true) {
  y_hatss <- purrr::map(
    sample_fn_list,
    \(sample_fn) mcmc_population(population, sample_fn, n_iter_true)
  )
  v_trues <- purrr::map_dbl(y_hatss, var)

  v_trues
}

compute_v_approxs <- function(population, v_approx_fn_list, x_names) {
  purrr::map_dbl(v_approx_fn_list, \(v_approx_fn) v_approx_fn(population, x_names))
}
