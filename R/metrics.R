#' Compute relative bias of estimations
#'
#' This computation requires producing multiple estimations and knowing the
#' true value (eventually by a Monte-Carlo approach).
#'
#' @param v_hats A vector containing multiple variance estimators.
#' @param v_true True variance.
#'
#' @return A scalar corresponding to the relative bias.
#' @export
#' @examples
#' compute_relative_bias(c(1,1,2,3,2), 2)
compute_relative_bias <- function(v_hats, v_true) {
  (mean(v_hats) - v_true) / v_true
}

#' Compute relative stability of estimations
#'
#' This computation requires producing multiple estimations and knowing the
#' true value (eventually by a Monte-Carlo approach).
#'
#' @param v_hats A vector containing multiple variance estimators.
#' @param v_true True variance.
#'
#' @return A scalar corresponding to the relative stability
#' @export
#' @examples
#' compute_relative_stability(c(1,1,2,3,2), 2)
compute_relative_stability <- function(v_hats, v_true) {
  sqrt(mean((v_hats - v_true)^2)) / v_true
}
