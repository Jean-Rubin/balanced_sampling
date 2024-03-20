lambdas <- list(
  c(1, -1),
  c(1, 2, -3),
  c(1, -1, 3),
  c(0, 1, -1),
  c(0, -1, 0)
)

describe("init_p()", {
  eps <- 1e-11
  for (lambda in lambdas) {
    p <- init_p(lambda, eps)

    it("gives a vector of the same size as lambda", {
      expect_equal(length(p), length(lambda))
    })
    it("gives a probability vector", {
      expect_equal(sum(p), 1)
      expect_true(all(p >= -eps))
    })
    it(add_context("is orthogonal to lambda", lambda), {
      expect_equal(sum(p * lambda), 0)
    })
    it("gives an interior solution, with strictly positive values" |>
      add_context(lambda),
      expect_true(all(p > 0) | all(lambda >= -eps) | all(lambda <= eps))
    )
  }
})

test_that("init_p() gives an interior solution", {
  eps <- 1e-11
  for (lambda in lambdas) {
    p <- init_p(lambda, eps)

    expect_equal(length(p), length(lambda))
    expect_equal(sum(p), 1)
    expect_equal(sum(p * lambda), 0)
    expect_true(
      all(lambda <= eps) | all(lambda >= eps) |
        all(p > 0 | lambda == 0)
    )
  }
})

test_that("max_entropy_p() gives a valid probability with a high entropy", {
  for (lambda in lambdas) {
    p0 <- init_p(lambda)
    p <- max_entropy_p(lambda)

    expect_equal(length(p), length(lambda))
    expect_equal(sum(p), 1)
    expect_equal(sum(p * lambda), 0)
    expect_gte(
      -sum(p * log(p), na.rm = TRUE),
      -sum(p0 * log(p0), na.rm = TRUE)
    )
  }
})

