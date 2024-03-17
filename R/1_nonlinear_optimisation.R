#' Produce an initial probability vector for the max-entropy sampling
#'
#' Every probability is initialized with a non-zero value.
#'
#' @param lambda
#' @param eps Threshold to consider a value as zero.
#'
#' @return A vector of probability with the same dimension as lambda.
#' @export
#'
#' @examples
#' init_p(c(1, 2, -3))
init_p <- function(lambda, eps = 1e-11) {
  p <- rep(0, length(lambda))
  lambda_neg <- which(lambda < -eps)
  lambda_pos <- which(lambda > eps)
  lambda_zero <- which(abs(lambda) <= eps)

  p[lambda_neg] <- sum(lambda[lambda_pos])
  p[lambda_pos] <- -sum(lambda[lambda_neg])
  p[lambda_zero] <- 1L

  p / sum(p)
}

#' Produce the max-entropy probabilities associated to a jumping step
#'
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
#' max_entropy_p(c(1, 2, -3), trace = 1)
max_entropy_p <- function(lambda, trace = 0) {
  p0 <- init_p(lambda)
  res <- Rsolnp::solnp(
    pars = p0,
    fun = \(p) sum(p * log(p), na.rm = TRUE),
    eqfun = \(p) c(sum(lambda * p), sum(p) - 1),
    eqB = c(0, 0),
    LB = rep(0, length(p0)),
    control = list(trace = trace)
  )

  res$pars
}

#' Sample an index with max-entropy under martingale constraint
#'
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
#' sample_max_entropy(c(1, 2, -3))
sample_max_entropy <- function(lambda) {
  p <- max_entropy_p(lambda)

  sample.int(length(lambda), size = 1, prob = p)
}
