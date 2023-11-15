population_test <- tibble::tibble(
  x = c(2, 2, 3, 5, 6, 7),
  y = c(5, 5, 6, 13, 12, 15),
  pi_i = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  pi_i_aux = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
)
sample_test_size <- sum(population_test$pi_i)


test_that("Balancing equation is valid after flight phase", {
  eps <- 1e-11
  X <- cbind(population_test$x, population_test$pi_i_aux)
  nb_aux_var <- dim(X)[2]
  s <- flight_wr(X, population_test$pi_i, eps)


  expect_equal(sum(abs(round(s) - s) > eps), nb_aux_var)
  expect_equal(sum(s), sample_test_size)
  expect_equal(sum(population_test$x * s / population_test$pi_i), sum(population_test$x))
})
