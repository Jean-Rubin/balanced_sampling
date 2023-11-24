#' Constraint a vector to fit inside bounds
#'
#' @param x A vector of numerical values.
#' @param thresh_min Lower bound.
#' @param thresh_max Upper bound.
#'
#' @return The same vector but whose values are clamped to the given bounds.
#' @export
#'
#' @examples clamp(c(0.1, 0.4, 0.9), thresh_min = 0.3, thresh_max = 0.8)
clamp <- function(x, thresh_min = 0.01, thresh_max = 0.99) {
  pmin.int(pmax.int(x, thresh_min), thresh_max)
}
