test_that("Clamp works on single values", {
  expect_equal(clamp(0.5, 0.1, 0.9), 0.5)
  expect_equal(clamp(0.5, 0.6, 0.9), 0.6)
  expect_equal(clamp(0.5, 0.1, 0.4), 0.4)
})

test_that("Clamp works on vector values", {
  expect_equal(clamp(c(0.1, 0.4, 0.9), 0.3, 0.8), c(0.3, 0.4, 0.8))
})
