#' Computation of the pseudo-inverse of a matrix
#'
#' @param x Matrix to pseudo-invert.
#' @param tol Threshold to consider an eigenvalue as null.
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' pseudo_inv(diag(c(1,2,3,0)), tol = 1e-3)
pseudo_inv <- function(x, tol = 1e-6) {
  svd_x <- svd(x)
  d <- diag(ifelse(svd_x$d < tol, 0, 1 / svd_x$d), nrow = length(svd_x$d))

  svd_x$v %*% d %*% t(svd_x$u)
}

#' Recover singular vectors corresponding to non-zero singular values
#'
#' @param X Matrix to recover singular values.
#' @param eps Threshold to consider a singular value as null.
#'
#' @return A reduced matrix.
#' @examples
#' reduc(matrix(c(1, 2, 3, 4), nrow = 2), eps = 1e-10)
reduc <- function(X, eps = 1e-11) {
  s <- svd(X)

  s$u[, s$d > eps, drop = FALSE]
}

#' Constraint a vector to fit inside bounds
#'
#' @param x A vector of numerical values.
#' @param thresh_min Lower bound.
#' @param thresh_max Upper bound.
#'
#' @return The same vector but whose values are clamped to the given bounds.
#' @export
#'
#' @examples
#' clamp(c(0.1, 0.4, 0.9), thresh_min = 0.3, thresh_max = 0.8)
clamp <- function(x, thresh_min = 0.01, thresh_max = 0.99) {
  pmin.int(pmax.int(x, thresh_min), thresh_max)
}