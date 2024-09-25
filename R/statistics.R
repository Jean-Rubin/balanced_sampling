#' Compute total of variable `y` of the population
#'
#' @param population A data frame corresponding to the population.
#'
#' @return The total.
#' @export
#' @examples
#' create_population(10) |> compute_total()
compute_total <- function(population) {
  population |>
    dplyr::summarize(y_true = sum(y)) |>
    dplyr::pull(y_true)
}

#' Covariance matrix of totals from multinomial sampling
#' 
#' @param x A matrix (n x p) representing a vector valued interest variable,
#'   where lines correspond to individual.
#' @param y A Matrix (n x q) representing a vector valued interest variable,
#'   where lines correspond to individual.
#' @param pi_i A vector (n) of inclusion probabilities of each individual.
#'
#' @return The associated covariance matrix (p x q).
#'
#' @examples
#' x <- matrix(c(1, 2, 1, 2), nrow = 2)
#' y <- matrix(c(1, 3, 4, 5, 6, 7), nrow = 2)
#' pi_i <- c(0.5, 0.8)
#' compute_covariance_multinomial(x, y, pi_i)
compute_covariance_multinomial <- function(x, y, pi_i) {
  n <- sum(pi_i)
  t_x <- colSums(x)
  t_y <- colSums(y)
  x_c <- sweep(x / pi_i, 2, t_x / n)
  y_c <- sweep(y / pi_i, 2, t_y / n)

  t(pi_i * x_c) %*% y_c
}

#' Multinomial approximation of variance of a balanced sample
#'
#' @param population A data frame corresponding to the population.
#'   It has:
#'   - `y`: the variable of interest.
#'   - `pi_i`: the inclusion probability of the individual.
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#'
#' @inheritParams pseudo_inv
#'
#' @return A variance approximation of the total estimator of `y`.
#' @export
compute_v_approx_multinomial <- function(population, x_names, tol = 1e-6) {
  x <- as.matrix(population[x_names])
  y <- as.matrix(population$y)
  v_x <- compute_covariance_multinomial(x, x, population$pi_i)
  c_xy <- compute_covariance_multinomial(x, y, population$pi_i)

  # pseudo inverse computation
  pseudo_inv_x <- pseudo_inv(v_x, tol)

  beta <- t(pseudo_inv_x) %*% c_xy
  e <- y - x %*% beta

  compute_covariance_multinomial(e, e, population$pi_i)
}

#' Deville-Tillé approximation of the variance of a balanced sample
#'
#' @param population A data frame corresponding to the population.
#'   It has:
#'   - `y`: the variable of interest.
#'   - `pi_i`: the inclusion probability of the individual.
#' @param x_names  A vector of the names of the auxiliary variables used to
#'   balance the sample.
#'
#' @inheritParams pseudo_inv
#'
#' @return A variance approximation of the total estimator of `y`.
#' @export
compute_v_approx_deville_tille <- function(population, x_names, tol = 1e-6) {
  x <- as.matrix(population[x_names])
  y <- as.matrix(population$y)
  pi_i <- population$pi_i
  n <- nrow(population)
  p <- length(x_names)

  x_diag_p <- t(x) %*% diag((1 - pi_i) / pi_i, ncol = length(pi_i))
  beta <- pseudo_inv(x_diag_p %*% x, tol) %*% x_diag_p %*% y
  e <- y - x %*% beta
  delta <- diag(pi_i * (1 - pi_i))

  (n / (n - p)) * (t(e / pi_i) %*% delta) %*% (e / pi_i)
}
