create_population <- function(n_tot, noise = 1) {
  tibble(
    const = 1,
    x1 = rnorm(n_tot, mean = 5, sd = 2),
    x2 = rnorm(n_tot, mean = -1, sd = 3),
    x3 = rnorm(n_tot, mean = 8, sd = 1),
    y = 2 * x1 + x2 - x3 + rnorm(n_tot, mean = 0, sd = noise)
  )
}

set_inclusion_proba <- function(population, pi_gen) {
  population |>
    mutate(
      pi_i = pi_gen(n()),
      pi_i_aux = pi_i
    )
}

pi_gen_beta <- function(a = 1, b = 5) {
  \(n) rbeta(n, a, b)
}

pi_gen_unif <- function(n_sample) {
  \(n) rep(n_sample / n, n)
}
