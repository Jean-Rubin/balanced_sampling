#' Repeated estimations of totals from independent samples
#'
#' This function generates multiple estimators of total obtained from repeated
#' sampling using the same procedure.
#'
#' @param population A data frame corresponding to the population.
#' @param sample_fn A sampling procedure as a function to apply to
#'   the population.
#'
#' @return A vector of estimated totals for each generated sample.
#' @export
#' @examples
#' population <- data.frame(y = 1:4, x1 = 5:8, pi_i = rep(0.2, 4))
#' mc_estimate_total(population, sampler_gen_base("x1"), 20)
mc_estimate_total <- function(
  population,
  sample_fn,
  n_iter,
  .progress = TRUE
) {
  purrr::map_dbl(
    seq_len(n_iter),
    \(x) {
      population |>
        sample_fn() |>
        estimate_total()
    },
    .progress = .progress
  )
}

#' Monte Carlo variance computation of estimators of total
#'
#' The variance associated to multiple sampling procedure can be applied.
#'
#' @param population A data frame corresponding to the population.
#' @param sample_fn_list A named list of sampling functions to apply to the
#'   population.
#' @param n_iter_true Number of iterations of the Monte Carlo simulation.
#'
#' @return A named vector with the variance associated to each sampling
#'   procedure. The names correspond to the names in the list.
#' @export
#' @examples
#' population <- data.frame(y = 1:4, x1 = 5:8, pi_i = rep(0.2, 4))
#' compute_true(population, list(base = sampler_gen_base("x1")), 20)
compute_v_trues <- function(population, sample_fn_list, n_iter_true) {
  y_hatss <- purrr::map(
    sample_fn_list,
    \(sample_fn) {
      mc_estimate_total(population, sample_fn, n_iter_true, .progress = FALSE)
    }
  )
  v_trues <- purrr::map_dbl(y_hatss, var)

  v_trues
}

#' Variance approximation of estimators of total
#'
#' Multiple variance approximations can be applied on the population.
#'
#' @param population A data frame corresponding to the population.
#' @param v_approx_fn_list A named list of sampling approximations to apply
#'   to the population.
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#'
#' @return A named vector with the variance approximations. The names
#'   correspond to the names in the list.
#' @export
#' @examples
#' population <- data.frame(y = 1:4, x1 = 5:8, pi_i = rep(0.2, 4))
#' compute_v_approxs(
#'   population,
#'   list(v_dt = compute_v_approx_deville_tille),
#'   "x1"
#' )
compute_v_approxs <- function(population, v_approx_fn_list, x_names) {
  purrr::map_dbl(
    v_approx_fn_list,
    \(v_approx_fn) v_approx_fn(population, x_names)
  )
}
