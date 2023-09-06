mcmc_population <- function(population, n_iter, .progress = TRUE) {
  purrr::map_dbl(
    seq_len(n_iter),
    \(x) population |>
      sample_population() |>
      estimate_total(),
    .progress = .progress
  )
}

