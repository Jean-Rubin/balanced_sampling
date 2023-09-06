create_population <- function(n_tot) {
  tibble(
    x = rnorm(n_tot, mean = 5, sd = 2),
    y = 2 * x + rnorm(n_tot, mean = 0, sd = 1)
  ) |>
    mutate(
      pi_i = rbeta(n(), 1, 5)
    )
}

