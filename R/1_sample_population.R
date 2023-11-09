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
        s_i = sampling::srswor(n_sample, n())
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
        s_i = sampling::samplecube(x, pi_i, comment = FALSE)
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
    wr_population <- tidyr::uncount(population, n_copy) |>
      mutate(
        pi_i = pi_i / n_copy,
        y = y / n_copy
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
        s_i = sample_cube(x, population[["pi_i"]])
      )
  }
}


sampler_gen_fuzzy_wr <- function(x_names, ...) {
  function(population) {
    x_local <- as.matrix(population[x_names])
    pi_i_local <- population[["pi_i"]]
    population |>
      mutate(
        s_i = sample_cube_wr(x_local, pi_i_local)
      )
  }
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

sample_cube_wr <- function(X, pik, eps = 1e-11) {
  algofastflightcube <- function(X, pik) {
    jump <- function(X, pik) {
      N = length(pik)
      p = round(length(X) / length(pik))
      X1 = cbind(X, rep(0, times = N))
      kern <- svd(X1)$u[, p + 1]
      listek = abs(kern) > eps
      # buff1 <- (1 - pik[listek]) / kern[listek]
      buff2 <- -pik[listek] / kern[listek]
      la1 <- min(buff2[buff2 > 0])
      pik1 <- pik + la1 * kern
      # buff1 <- -(1 - pik[listek]) / kern[listek]
      buff2 <- pik[listek] / kern[listek]
      la2 <- min(buff2[buff2 > 0])
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
    n_select <- numeric(length = N)
    # X <- array(X, c(N, p))
    A <- X / pik
    psik <- pik[1:(p + 1)]
    B <- A[1:(p + 1), ]
    ind <- 1:(p + 1)
    # B <- array(B, c(p + 1, p))
    pp = p + 2
    liste_rej <- 1:(p + 1)
    liste_next <- 1:(p + 1)
    while (pp <= N + 1) {
      psik[liste_rej] <- pik[liste_next]
      B[liste_rej, ] <- A[liste_next, ]
      ind[liste_rej] <- liste_next

      psik <- abs(jump(B, psik))
      # liste_1 <- psik > (1 - eps)
      # liste_0 <- psik < eps
      # liste <- liste_0 | liste_1
      n_select[ind] <- n_select[ind] + floor(psik)
      psik <- psik - floor(psik)
      liste_rej <- which(psik < eps)
      liste_next <- pp - 1 + seq_along(liste_rej)
      pp <- pp + length(liste_rej)

      # i <- 0
      # while (i <= p && pp <= N + 1) {
      #   i = i + 1
      #   if (liste[i]) {
      #     # pik[ind[i]] = psik[i]
      #     psik[i] = pik[pp]
      #     B[i, ] = A[pp, ]
      #     B = array(B, c(p + 1, p))
      #     ind[i] = pp
      #     pp = pp + 1
      #   }
      # }
    }
    if (length(pik[abs(pik - round(pik)) > eps]) == (p + 1))
      psik <- jump(B, psik)
    pik[ind] = psik

    res <- n_select
    res[ind] <- res[ind] + psik

    res
  }
  reduc <- function(X) {
    eps = 1e-11
    N = dim(X)[1]
    Re = svd(X)
    array(Re$u[, (Re$d > eps)] , c(N, sum(as.integer(Re$d > eps))))
  }

  N = length(pik)

  p = round(length(X) / length(pik))
  X <- array(X, c(N, p))
  o <- sample.int(N, N)
  liste <- o[abs(pik[o] - round(pik[o])) > eps]

  pikbon <- pik[liste]

  Nbon = length(pikbon)

  Xbon <- array(X[liste, ] , c(Nbon, p))
  pikstar <- pik
  if (Nbon > p) {
    pikstarbon <- algofastflightcube(Xbon, pikbon)
    pikstar[liste] = pikstarbon
  }

  pbon = dim(Xbon)[2]
  while (Nbon > pbon && Nbon > 0) {
    pikstarbon <- algofastflightcube(Xbon / pik[liste] * pikbon, pikbon)
    pikstar[liste] = pikstarbon
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
