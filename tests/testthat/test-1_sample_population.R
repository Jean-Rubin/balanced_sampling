eps <- 1e-11
population_test <- tibble::tibble(
  x = c(2, 2, 3, 5),
  y = c(3, 5, 6, 13),
  pi_i = c(0.5, 0.5, 0.5, 0.5),
  pi_i_aux = c(0.5, 0.5, 0.5, 0.5)
)
X <- cbind(population_test$x, population_test$pi_i_aux)
sample_test_size <- sum(population_test$pi_i)

test_that("Balancing equation is valid after flight phase", {
  nb_aux_var <- dim(X)[2]
  for (k in seq_len(100)) {
    s <- flight_wr(X, population_test$pi_i, eps)
    expect_equal(sum(abs(round(s) - s) > eps), nb_aux_var)
    expect_equal(sum(s), sample_test_size)
    expect_equal(sum(population_test$x * s / population_test$pi_i), sum(population_test$x))
  }
})

test_that("Balancing equation is valid after flight phase 2", {
  nb_aux_var <- dim(X)[2]
  s <- sample_cube_wr(X, population_test$pi_i, eps)

  expect_equal(sum(abs(round(s) - s) > eps), nb_aux_var)
  expect_equal(sum(s), sample_test_size)
  expect_equal(sum(population_test$x * s / population_test$pi_i), sum(population_test$x))
})



test_that("jump_wr_ent() satisfies basic properties", {
  for (k in seq_len(100)) {
    jump_result <- jump_wr_ent(X, population_test$pi_i, eps)
    n_select <- jump_result$pi_jump[jump_result$k_select]

    # The jumping coordinate is fixed at an integer
    expect_true(abs(n_select - round(n_select)) < eps)
    # Jumping probability is positive
    expect_true(all(jump_result$pi_jump > -eps))
    # Balancing equation is valid
    expect_equal(sum(population_test$x * jump_result$pi_jump / population_test$pi_i), sum(population_test$x))
  }
})

test_that("jump_wr_ent() works on edge case", {
  X <- matrix(c(1, 1, 1), ncol = 1)
  pik <- c(1e-13, 0.9, 0)
  pik_2 <- c(-1e-13, 0.9, 0)

  expect_no_error(jump_wr_ent(X, pik, eps = 1e-11))
  expect_no_error(jump_wr_ent(X, pik_2, eps = 1e-11))
})

test_that("flight_phase_wr_ent() satisfies basic properties", {
  s <- flight_phase_wr_ent(X, population_test$pi_i, eps)

  expect_equal(sum(s), sample_test_size)
  expect_equal(sum(population_test$x * s / population_test$pi_i), sum(population_test$x))
})
