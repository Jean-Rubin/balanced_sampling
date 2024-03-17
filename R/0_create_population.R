create_population <- function(n_tot, noise = 1) {
  tibble(
    const = 1,
    x1 = rnorm(n_tot, mean = 5, sd = 2),
    x2 = rnorm(n_tot, mean = -1, sd = 3),
    x3 = rnorm(n_tot, mean = 8, sd = 1),
    x4 = rnorm(n_tot, mean = 0, sd = noise),
    y = 2 * x1 + x2 - x3 + x4
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

# Note: this does not depend on n, it ideally should depend on population
# Now we will assume that the population is not dynamically changing, so it is
# possible to precompute the sizes.
pi_gen_prop_size <- function(n_sample, x) {
  \(n) sampling::inclusionprobabilities(x, n_sample)
}
