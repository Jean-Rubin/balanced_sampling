#' Standard without-replacement sampling
#'
#' @param population
#'
#' @return A sample with the same columns as the population, but with an
#'   indicator of selection `s_i`.
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
#' @param population
#'
#' @return A sample with the same columns as the population, but with an
#'   indicator of selection `s_i`.
#' @export
sampler_gen <- function(x_names, ...) {
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
#' @param population
#'
#' @return A sample with the same columns as the population, but with an
#'   indicator of selection `s_i`.
#' @export
sampler_gen_wr_copy <- function(x_names, ...) {
  function(population) {
    n_copy <- sum(population[["pi_i"]])
    wr_population <- population |>
      tidyr::uncount(.env$n_copy) |>
      mutate(
        pi_i = pi_i / .env$n_copy,
        y = y / .env$n_copy
      )

    sampler_gen(x_names)(wr_population)
  }
}

#' Cube method fuzzy sampling (flight phase)
#'
#' @param population
#'
#' @return A sample with the same columns as the population, but with an
#'   indicator of selection `s_i` that can be between 0 and 1.
#' @export
sampler_gen_fuzzy_custom <- function(x_names,  ...) {
  function(population) {
    x <- as.matrix(population[x_names])
    population |>
      mutate(
        s_i = sample_cube(.env$x, pi_i)
      )
  }
}

sampler_gen_fuzzy_wr <- function(x_names, ...) {
  function(population) {
    x <- as.matrix(population[x_names])
    population |>
      mutate(
        s_i = sample_cube_wr(.env$x, pi_i)
      )
  }
}

reduc <- function(X, eps = 1e-11) {
  N <- dim(X)[1]
  s <- svd(X)

  s$u[, s$d > eps, drop = FALSE]
}


sample_cube <- function(X, pik, eps = 1e-11) {
  algofastflightcube <- function(X, pik) {
    jump <- function(X, pik) {
      N = length(pik)
      p = round(length(X) / length(pik))
      # X <- array(X, c(N, p))
      X1 = cbind(X, rep(0, times = N))
      kern <- svd(X1)$u[, p + 1]
      listek = abs(kern) > eps
      buff1 <- (1 - pik[listek]) / kern[listek]
      buff2 <- -pik[listek] / kern[listek]
      la1 <- min(c(buff1[(buff1 > 0)] , buff2[(buff2 > 0)]))
      pik1 <- pik + la1 * kern
      buff1 <- -(1 - pik[listek]) / kern[listek]
      buff2 <- pik[listek] / kern[listek]
      la2 <- min(c(buff1[(buff1 > 0)] , buff2[(buff2 > 0)]))
      pik2 <- pik - la2 * kern
      q <- la2 / (la1 + la2)
      if (runif(1) < q)
        pikn <- pik1
      else
        pikn <- pik2
      pikn
    }

    N = length(pik)
    p = round(length(X) / length(pik))
    # X <- array(X, c(N, p))
    A <- X / pik
    B <- A[1:(p + 1), ]
    psik <- pik[1:(p + 1)]
    ind <- seq(1, p + 1, 1)
    pp = p + 2
    B <- array(B, c(p + 1, p))
    while (pp <= N) {
      psik <- jump(B, psik)
      liste <- (psik > (1 - eps) | psik < eps)
      i <- 0
      while (i <= (p) && pp <= N) {
        i = i + 1
        if (liste[i]) {
          pik[ind[i]] = psik[i]
          psik[i] = pik[pp]
          B[i, ] = A[pp, ]
          B = array(B, c(p + 1, p))
          ind[i] = pp
          pp = pp + 1
        }
      }
    }
    if (length(pik[(pik > eps & pik < (1 - eps))]) == (p + 1))
      psik <- jump(B, psik)
    pik[ind] = psik
    pik
  }
  reduc <- function(X) {
    eps = 1e-11
    N = dim(X)[1]
    Re = svd(X)
    array(Re$u[, (Re$d > eps)] , c(N, sum(as.integer(Re$d > eps))))
  }

  N = length(pik)

  p = round(length(X) / length(pik))
  # X <- array(X, c(N, p))
  o <- sample.int(N, N)
  liste <- o[(pik[o] > eps & pik[o] < (1 - eps))]

  pikbon <- pik[liste]

  Nbon = length(pikbon)

  Xbon <- array(X[liste, ] , c(Nbon, p))
  pikstar <- pik
  if (Nbon > p) {
    pikstarbon <- algofastflightcube(Xbon, pikbon)
    pikstar[liste] = pikstarbon
  }
  # liste <- o[(pikstar[o] > eps & pikstar[o] < (1 - eps))]
  # pikbon <- pikstar[liste]
  # Nbon = length(pikbon)
  # Xbon <- array(X[liste, ] , c(Nbon, p))
  # pbon = dim(Xbon)[2]
  # if (Nbon > 0) {
  #   Xbon = reduc(Xbon)
  #   pbon = dim(Xbon)[2]
  # }

  pbon = dim(Xbon)[2]
  while (Nbon > pbon && Nbon > 0) {
    pikstarbon <- algofastflightcube(Xbon / pik[liste] * pikbon, pikbon)
    pikstar[liste] = pikstarbon
    liste <- o[(pikstar[o] > eps & pikstar[o] < (1 - eps))]
    pikbon <- pikstar[liste]
    Nbon = length(pikbon)
    Xbon <- array(X[liste, ] , c(Nbon, p))
    if (Nbon > 0) {
      Xbon = reduc(Xbon)
      pbon = dim(Xbon)[2]
    }
  }

  pikstar
}

jump_wr <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  X1 <- cbind(X, rep(0, times = N))
  kern <- svd(X1)$u[, p + 1]
  idx_kern <- abs(kern) > eps
  buff <- -pik[idx_kern] / kern[idx_kern]
  la1 <- min(buff[buff > 0])
  la2 <- min(-buff[buff < 0])
  q <- la1 / (la1 + la2)

  pik + (la1 - (la1 + la2) * rbinom(1, 1, q)) * kern
}

flight_wr <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  A <- X / pik
  n_select <- numeric(N) # list of complete selection of units
  psik <- numeric(p + 1) # buffer of probability of selecting a unit
  B <- matrix(0, nrow = p + 1, ncol = p)
  ind <- 1:(p + 1)
  idx_rej <- 1:(p + 1)
  idx_next <- 1:(p + 1)
  pp <- p + 1 # last p
  while (pp <= N) {
    psik[idx_rej] <- pik[idx_next]
    B[idx_rej, ] <- A[idx_next, ]
    ind[idx_rej] <- idx_next

    # Absolute value to prevent numerical instability of sign near zero
    # This matter when calling the floor function
    psik <- abs(jump_wr(B, psik, eps))

    n_select[ind] <- n_select[ind] + floor(psik)
    psik <- psik - floor(psik)
    idx_rej <- which(psik < eps)
    idx_next <- pp + seq_along(idx_rej)
    pp <- pp + length(idx_rej)
  }
  if (length(pik[abs(pik - round(pik)) > eps]) == (p + 1))
    psik <- jump_wr(B, psik, eps)

  # We regroup selected units and their probabilistic part
  n_select[ind] <- n_select[ind] + psik

  n_select
}

sample_cube_wr <- function(X, pik, eps = 1e-11) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  o <- sample.int(N, N)
  liste <- o[abs(pik[o] - round(pik[o])) > eps]

  pikbon <- pik[liste]

  Nbon <- length(pikbon)
  Xbon <- X[liste, ]
  pikstar <- pik
  if (Nbon > p) {
    pikstarbon <- algofastflightcube(Xbon, pikbon)
    pikstar[liste] <- pikstarbon
  }

  pbon <- dim(Xbon)[2]
  while (Nbon > pbon && Nbon > 0) {
    pikstarbon <- algofastflightcube(Xbon / pik[liste] * pikbon, pikbon)
    pikstar[liste] <- pikstarbon
    liste <- o[abs(pikstar[o] - round(pikstar[o])) > eps]
    pikbon <- pikstar[liste]
    Nbon = length(pikbon)
    Xbon <- array(X[liste, ] , c(Nbon, p))
    if (Nbon > 0) {
      Xbon = reduc(Xbon)
      pbon = dim(Xbon)[2]
    }
  }

  pikstar
}
