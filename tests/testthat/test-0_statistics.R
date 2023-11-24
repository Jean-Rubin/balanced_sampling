population_test <- tibble::tibble(
  x = c(2, 2, 3, 5),
  y = c(4, 2, 4, 6),
  pi_i = c(0.5, 0.5, 0.5, 0.5),
  pi_i_aux = c(0.5, 0.5, 0.5, 0.5)
)

test_that("Computation of total works", {
  expect_equal(compute_total(population_test), 16)
})

test_that("Pseudo inverse return a matrix", {
  expect_equal(pseudo_inv(matrix(2)), matrix(0.5))
})

test_that("Multinomial covariance computation works", {
  x <- as.matrix(population_test["x"])
  y <- as.matrix(population_test["y"])
  pi_i <- population_test$pi_i
  expect_equal(
    compute_covariance_multinomial(x, x, pi_i),
    matrix(12, dimnames = list("x", "x"))
  )
  expect_equal(
    compute_covariance_multinomial(y, y, pi_i),
    matrix(16, dimnames = list("y", "y"))
  )
  expect_equal(
    compute_covariance_multinomial(x, y, pi_i),
    matrix(12, dimnames = list("x", "y"))
  )
})

test_that("Multinomial approximation works", {
  expect_equal(
    compute_v_approx_multinomial(population_test, "x"),
    matrix(4)
  )
  expect_equal(
    compute_v_approx_multinomial(population_test, "pi_i_aux"),
    matrix(16)
  )
  expect_equal(
    compute_v_approx_multinomial(population_test, c("x", "pi_i_aux")),
    matrix(4)
  )
})
