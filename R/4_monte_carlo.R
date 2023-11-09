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

# compute_var_stat_n <- function(population, n_s) {
#   population <- population |> set_inclusion_proba(pi_gen_unif(n_s))
#   compute_v_trues(population, sample_fn_list, mcmc_population, n_iter_true)
# }
#
# compute_var_stat_fuzzy <- function(n_s) {
#   population <- population |> set_inclusion_proba(pi_gen_unif(n_s))
#   y_hats <- mcmc_fuzzy(population, sample_cube, n_iter_true)
#   y_hats_wr <- mcmc_fuzzy(population, sample_cube_wr, n_iter_true)
#
#   v_hat <- var(y_hats)
#   v_hat_wr <- var(y_hats_wr)
#
#   list(
#     v_hat = v_hat,
#     v_hat_wr = v_hat_wr
#   )
# }
