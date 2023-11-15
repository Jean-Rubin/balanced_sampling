population_test <- tibble::tibble(
  x = c(2, 2, 3, 5, 6, 7),
  y = c(5, 5, 6, 13, 12, 15),
  pi_i = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  pi_i_aux = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
)

test_that("Computation of total works", {
  expect_equal(compute_total(population_test), 56)
})
