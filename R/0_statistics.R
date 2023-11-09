#' Calcul du pseudo-inverse d'une matrice
#'
#' @param x Matrice à pseudo-inverser.
#' @param tol Seuil de tolérance pour déterminer les valeurs propres nulles.
#'
#' @return Une matrice.
#' @export
#'
#' @examples pseudo_inv(diag(c(1,2,3,0)), tol = 1e-3)
pseudo_inv <- function(x, tol = 1e-6) {
  svd_x <- svd(x)
  d <- diag(ifelse(svd_x$d < tol, 0, 1 / svd_x$d))

  svd_x$v %*% d %*% t(svd_x$u)
}

#' Calcul du total sur la variable `y` de la population
#'
#' @param population Population
#'
#' @return Le total
#' @export
compute_total <- function(population) {
  population |>
    summarize(y_true = sum(y)) |>
    pull(y_true)
}

#' Calcul de la matrice de variance-covariance de deux totaux de variables
#' d'intérêt associées à un tirage multinomial
#'
#' @param x Matrice (n x p) représentant une variable d'intérêt vectorielle, les
#'   lignes correspondant aux individus.
#' @param y Matrice (n x q) représentant une variable d'intérêt vectorielle, les
#'   lignes correspondant aux individus.
#' @param pi_i Vecteur des probabilités d'inclusion de chaque individu.
#'
#' @return Une matrice (p x q) de variance-covariance
compute_covariance_multinomial <- function(x, y, pi_i) {
  n <- sum(pi_i)
  t_x <- colSums(x)
  t_y <- colSums(y)
  x_c <- sweep(x / pi_i, 2, t_x / n)
  y_c <- sweep(y / pi_i, 2, t_y / n)

  t(pi_i * x_c) %*% y_c
}

#' Approximation multinomiale de la variance d'un tirage équilibré
#'
#' @param population Population.
#' @param x_names Vecteur des noms des variables auxiliaires utilisées pour
#'   l'équilibrage.
#' @inheritParams pseudo_inv
#'
#' @return Une approximation de la variance du total sur la variable `y` de la
#'   population.
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

#' Approximation de Deville-Tillé de la variance d'un tirage équilibré
#'
#' @param population Population.
#' @param x_names Vecteur des noms des variables auxiliaires utilisées pour
#'   l'équilibrage.
#' @inheritParams pseudo_inv
#'
#' @return Une approximation de la variance du total sur la variable `y` de la
#'   population.
#' @export
compute_v_approx_deville_tille <- function(population, x_names, tol = 1e-6) {
  x <- as.matrix(population[x_names])
  y <- as.matrix(population$y)
  pi_i <- population$pi_i
  n <- nrow(population)
  p <- length(x_names)

  x_diag_p <- t(x) %*% diag((1 - pi_i) / pi_i)
  beta <- pseudo_inv(x_diag_p %*% x, tol) %*% x_diag_p %*% y
  e <- y - x %*% beta
  delta <- diag(pi_i * (1 - pi_i))

  (n / (n - p)) * (t(e / pi_i) %*% delta) %*% (e / pi_i)
}
