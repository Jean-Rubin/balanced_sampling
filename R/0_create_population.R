# Population creation ---------------------------------------------------------

#' Create a population
#'
#' @param n_tot Total number of individual.
#' @param noise Noise parameter in the generation of the variable of interest.
#'
#' @return A `tibble` data frame.
#' @export
#' @examples
#' create_population(10, noise = 1)
create_population <- function(n_tot, noise = 1) {
  tibble(
    const = 1,
    x1 = rnorm(n_tot, mean = 10, sd = 1),
    x2 = rnorm(n_tot, mean = -1, sd = 5),
    x3 = rnorm(n_tot, mean = 3, sd = 1),
    eps = rnorm(n_tot, mean = 0, sd = noise),
    y = 10 + 4 * x1 + x2 - 3 * x3 + eps
  )
}

## Generator helper -----------------------------------------------------------

#' Set the inclusion probabilities of the population
#'
#' Description
#'
#' @param population A data frame corresponding to the population.
#' @param pi_gen A probability generator. This generator must take
#'   the number of individual as a parameter.
#'
#' @return The same data frame with two additionnal columns:
#'   `pi_i` and `pi_i_aux`.
#' @examples
#' create_population(10) |> set_inclusion_proba(pi_gen_unif)
set_inclusion_proba <- function(population, pi_gen) {
  population |>
    mutate(
      pi_i = pi_gen(n()),
      pi_i_aux = pi_i
    )
}

### Inclusion probability generator -------------------------------------------

#' Beta distribution probability generator
#'
#' @param a Parameter of the beta distribution.
#' @param b Parameter of the beta distribution.
#'
#' @return A probability generator.
#' @export
#' @examples
#' pi_gen_beta(1, 5)
pi_gen_beta <- function(a = 1, b = 5) {
  \(n) rbeta(n, a, b)
}

#' Uniform distribution probability generator
#'
#' @param n_sample Number of element in the sample.
#'
#' @return A probability generator.
#' @export
#' @examples
#' pi_gen_unif(10)
pi_gen_unif <- function(n_sample) {
  \(n) rep(n_sample / n, n)
}

#' Proportional to variable - probability generator
#'
#' Note: this does not depend on n, it ideally should depend on population
#' Now we will assume that the population is not dynamically changing, so it is
#' possible to precompute the auxiliary variables.
#'
#' @param n_sample Number of element in the sample.
#' @param x Auxiliary variable to which the probabilities should be
#'   proportional to.
#'
#' @return A probability generator.
#' @export
#' @examples
#' # This generator assumes that the population is of size 7
#' pi_gen_prop_size(2, c(1,2,3,5,2,1,1))
pi_gen_prop_size <- function(n_sample, x) {
  \(n) sampling::inclusionprobabilities(x, n_sample)
}
