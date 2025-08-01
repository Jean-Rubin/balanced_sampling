# Sampling strategies ---------------------------------------------------------

#' Standard without-replacement sampling
#'
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i`.
sampler_gen_srswor <- function(...) {
  function(population) {
    n_sample <- round(sum(population[["pi_i"]]))
    population |>
      mutate(
        s_i = sampling::srswor(.env$n_sample, n())
      )
  }
}

#' Standard balanced sampling using the cube method
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i`.
#' @export
sampler_gen_base <- function(x_names, ...) {
  function(population) {
    x <- as.matrix(population[x_names])
    population |>
      mutate(
        s_i = sampling::samplecube(.env$x, pi_i, comment = FALSE)
      ) |>
      filter(s_i == 1) # Optimization by reducing data size
  }
}

#' With-replacement cube method sampling using dimensional increase
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param y_names A vector with the names of the variables of interest.
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i` that can be between 0 and 1.
#' @export
sampler_gen_wr_copy <- function(x_names, y_names = "y", ...) {
  function(population) {
    # Random rounding if non-integer sum
    n_copy <- floor(runif(1) + sum(population[["pi_i"]]))
    wr_population <- population |>
      tidyr::uncount(.env$n_copy) |>
      mutate(
        across(all_of(y_names), \(y) y / .env$n_copy),
        pi_i = pi_i / .env$n_copy
      )

    sampler_gen_base(x_names)(wr_population)
  }
}

#' With-replacement cube method sampling using dimensional increase (flight phase)
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param y_names A vector with the names of the variables of interest.
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i` that can be between 0 and 1.
#' @export
sampler_gen_flight_wr_copy <- function(x_names, y_names = "y", ...) {
  function(population) {
    # Random rounding if non-integer sum
    n_copy <- floor(runif(1) + sum(population[["pi_i"]]))
    wr_population <- population |>
      tidyr::uncount(.env$n_copy) |>
      mutate(
        across(all_of(y_names), \(y) y / .env$n_copy),
        pi_i = pi_i / .env$n_copy
      )

    sampler_gen_flight_base(x_names)(wr_population)
  }
}


#' Cube method (flight phase)
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i` that can be between 0 and 1.
#' @export
sampler_gen_flight_base <- function(x_names,  ...) {
  function(population) {
    x <- as.matrix(population[x_names])
    population |>
      mutate(
        s_i = sampling::fastflightcube(.env$x, pi_i, comment = FALSE)
      )
  }
}

#' With-replacement exhaustion cube method (flight phase)
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i` that can be between 0 and 1.
#' @export
sampler_gen_flight_wr_exh <- function(x_names, ...) {
  function(population) {
    x <- as.matrix(population[x_names])
    population |>
      mutate(
        s_i = flight_wr_exh(.env$x, pi_i)
      )
  }
}

#' With-replacement stepwise max-entropy cube method (flight phase)
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param ... Unused extra arguments.
#'
#' @return A sampling function, which given a population, returns a sample with
#'   the same columns as the population, but with an indicator of selection
#'   named `s_i` that can be between 0 and 1.
#' @export
sampler_gen_flight_wr_ent <- function(x_names, ...) {
  function(population) {
    x <- as.matrix(population[x_names])
    population |>
      mutate(
        s_i = flight_wr_ent(.env$x, pi_i)
      )
  }
}


## With-replacement exhaustion cube method ------------------------------------

#' Jump step in the with-replacement exhaustion cube method
#'
#' Will select a random direction in the kernel of the balancing matrix and move
#' the current sampling/inclusion probability vector to a new state. Given a
#' direction, two candidate states are selected, corresponding to the moment
#' where a unit is unselected.
#'
#' @param X The balancing matrix.
#' @param pik The sampling/inclusion probability vector before jumping.
#' @param eps A threshold for determining if a value is null.
#'
#' @return The sampling vector after the jump.
#'
#' @examples
jump_wr_exh <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  X1 <- cbind(X, rep(0, times = N))
  kern <- svd(X1)$u[, p + 1]
  idx_kern <- abs(kern) > eps
  buff <- -pik[idx_kern] / kern[idx_kern]
  la1 <- min(buff[buff >= 0])
  la2 <- min(-buff[buff <= 0])
  q <- la1 / (la1 + la2 + eps)

  pik + (la1 - (la1 + la2) * rbinom(1, 1, q)) * kern
}

#' Flight phase of the with-replacement exhaustion cube method
#'
#' @param X The balancing matrix.
#' @param pik The inclusion probability vector.
#' @param eps A threshold for determining if a value is null.
#'
#' @return A sampling vector of multiplicities
#' @export
#'
#' @examples
flight_wr_exh <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  random_perm <- sample.int(N, N)
  A <- (X / pik)[random_perm, ]
  n_select <- numeric(N) # list of complete selection of units
  psik <- numeric(p + 1) # buffer of probability of selecting a unit
  B <- matrix(0, nrow = p + 1, ncol = p)
  ind <- 1:(p + 1)
  idx_rej <- 1:(p + 1)
  idx_next <- 1:(p + 1)
  while (length(idx_next) > 0 && idx_next[1] <= N) {
    nb_valid_ind <- sum(idx_next <= N)
    idx_next <- idx_next[seq_len(nb_valid_ind)]
    idx_rej <- idx_rej[seq_len(nb_valid_ind)]
    psik[idx_rej] <- pik[idx_next]
    B[idx_rej, ] <- A[idx_next, ]
    ind[idx_rej] <- idx_next

    # Absolute value to prevent numerical instability of sign near zero
    # This matter when calling the floor function
    psik <- abs(jump_wr_exh(B, psik, eps))

    n_select[ind] <- n_select[ind] + floor(psik + eps)
    psik <- psik - floor(psik + eps)
    idx_rej <- which(psik < eps)
    idx_next <- idx_next[length(idx_next)] + seq_along(idx_rej)
  }

  # We regroup selected units and their probabilistic part
  n_select[ind] <- n_select[ind] + psik

  n_select[order(random_perm)]
}

sample_cube_wr_exh <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  o <- sample.int(N, N)

  ind_remain <- o[abs(pik[o] - round(pik[o])) > eps]
  X_remain <- X[ind_remain, , drop = FALSE]
  while (length(ind_remain) > dim(X_remain)[2]) {
    pik[ind_remain] <- flight_wr_exh(X_remain, pik[ind_remain], eps)
    ind_remain <- o[abs(pik[o] - round(pik[o])) > eps]
    X_remain <- reduc(X[ind_remain, , drop = FALSE])
  }

  pik
}


## With replacement cube with max-entropy jump --------------------------------

jump_wr_ent <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  X1 <- cbind(X, rep(0, times = N))
  u_kern <- svd(X1)$u[, p + 1]

  # find jump candidates
  # Some values are eventually infinite here
  lambdak_0 <- -pik / u_kern
  lambda_min <- max(lambdak_0[u_kern > eps]) - eps
  lambda_max <- min(lambdak_0[u_kern < -eps]) + eps

  lambda_low <- (lambda_max + (u_kern > 0) * (lambda_min - lambda_max))
  mk_min <- ceiling(pik + u_kern * lambda_low)
  mk_max <- floor(pik + u_kern * (lambda_min + lambda_max - lambda_low))
  mk <- purrr::map2(mk_min, mk_max, \(x, y) if (x > y) NULL else seq(x, y))

  # Compute the associated lambda and pick one
  lambdak <- purrr::pmap(
    list(mk, pik, u_kern),
    \(m, p, u) (m - p) / u
  ) |>  purrr::list_c()
  k_selects <- rep(seq_along(mk), purrr::map_dbl(mk, length))
  lambdak

  idx_select <- sample_max_entropy(lambdak)
  k_select <- k_selects[idx_select]
  lambda_select <- lambdak[idx_select]

  list(
    k_select = k_select,
    pi_jump = pik + lambda_select * u_kern
  )
}

flight_wr_ent <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  A <- X / pik

  # Final pi after flight phase
  #   will store number of times we select a unit
  n_select <- numeric(N)
  # Jumping pi, its length decreases as we select units
  pik_jump <- pik
  # Indices of units that are not selected yet
  idx_remain <- seq_len(N)
  while (length(idx_remain) > p) {
    B <- A[idx_remain, ]
    jump_result <- jump_wr_ent(B, pik_jump, eps)
    pik_jump <- jump_result$pi_jump
    k_select <- jump_result$k_select

    n_select[idx_remain][k_select] <- pik_jump[k_select]
    idx_remain <- idx_remain[-k_select]
    pik_jump <- pik_jump[-k_select]
  }

  # We regroup selected units and their probabilistic part
  n_select[idx_remain] <- n_select[idx_remain] + pik_jump

  n_select
}
