eps <- 1e-11
population_test <- tibble::tibble(
  x = c(2, 2, 3, 5),
  y = c(3, 5, 6, 13),
  pi_i = c(0.5, 0.5, 0.5, 0.5),
  pi_i_aux = c(0.5, 0.5, 0.5, 0.5)
)
X <- cbind(population_test$x, population_test$pi_i_aux)
sample_test_size <- sum(population_test$pi_i)

describe("flight_wr_exh()", {
  nb_aux_var <- dim(X)[2]
  set.seed(1234)
  for (k in seq_len(100)) {
    s <- flight_wr_exh(X, population_test$pi_i, eps)
    it("has at most p non-integer coordinates, where p is the number of auxiliary variables" |> add_context(k), {
      expect_lte(sum(abs(round(s) - s) > eps), nb_aux_var)
    })
    it("verifies the balancing constraints" |> add_context(k), {
      expect_equal(colSums(X * s / population_test$pi_i), colSums(X))
    })
  }
})

describe("sample_cube_wr_exh()", {
  nb_aux_var <- dim(X)[2]
  set.seed(1234)
  for (k in seq_len(100)) {
    s <- sample_cube_wr_exh(X, population_test$pi_i, eps)
    it("has at most p non-integer coordinates, where p is the number of auxiliary variables" |>
        add_context(k),
      expect_lte(sum(abs(round(s) - s) > eps), nb_aux_var)
    )
    it("verifies the balancing constraints" |>
        add_context(k),
      expect_equal(colSums(X * s / population_test$pi_i), colSums(X))
    )
  }
})

describe("jump_wr_ent()", {
  set.seed(1234)
  for (k in seq_len(100)) {
    jump_result <- jump_wr_ent(X, population_test$pi_i, eps)
    n_select <- jump_result$pi_jump[jump_result$k_select]

    it("fixes the jumping coordinate to an integer" |>
        add_context(k),
      expect_true(abs(n_select - round(n_select)) < eps)
    )
    it("gives positive jumping probabilities" |>
        add_context(k),
      expect_true(all(jump_result$pi_jump > -eps))
    )
    it("verifies the balancing constraints" |>
        add_context(k),
      expect_equal(colSums(X * jump_result$pi_jump / population_test$pi_i), colSums(X))
    )
  }

  it("works on edge case" |> add_context(k), {
    X <- matrix(c(1, 1, 1), ncol = 1)
    pik <- c(1e-13, 0.9, 0)
    pik_2 <- c(-1e-13, 0.9, 0)

    expect_no_error(jump_wr_ent(X, pik, eps = 1e-11))
    expect_no_error(jump_wr_ent(X, pik_2, eps = 1e-11))
  })
})

describe("flight_wr_ent()", {
  set.seed(1234)
  s <- flight_wr_ent(X, population_test$pi_i, eps)

  it("verifies the balancing constraints", {
    expect_equal(colSums(X * s / population_test$pi_i), colSums(X))
  })
})
