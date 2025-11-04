eps <- 1e-11
# Population test -------------------------------------------------------------
# Equal probability
population_test_1 <- tibble::tibble(
  x = c(2, 2, 3, 5),
  y = c(3, 5, 6, 13),
  pi_i = c(0.5, 0.5, 0.5, 0.5),
  pi_i_aux = c(0.5, 0.5, 0.5, 0.5)
)
X_1 <- cbind(population_test_1$x, population_test_1$pi_i_aux)

# Unequal probability, not an integer sample size
population_test_2 <- tibble::tibble(
  x = c(2, 2, 3, 5),
  y = c(3, 6, 6, 13),
  pi_i = c(0.2, 0.8, 0.1, 0.4),
  pi_i_aux = c(0.2, 0.8, 0.1, 0.4)
)
X_2 <- cbind(population_test_2$x, population_test_2$pi_i_aux)

# Tests -----------------------------------------------------------------------

describe("jump_wr_exh()", {
  set.seed(1234)
  for (k in seq_len(10)) {
    jump_result_1 <- jump_wr_exh(X_1 / population_test_1$pi_i, population_test_1$pi_i, eps)

    it("fixes the jumping coordinate to an integer" |> add_context(k),
      expect_gt(
        sum(abs(jump_result_1 - round(jump_result_1)) < eps),
        sum(abs(population_test_1$pi_i - round(population_test_1$pi_i)) < eps)
      )
    )
    it("gives positive jumping probabilities" |> add_context(k),
      expect_true(all(jump_result_1 > -eps))
    )
    it("verifies the balancing constraints" |> add_context(k),
      expect_equal(colSums(X_1 * jump_result_1 / population_test_1$pi_i), colSums(X_1))
    )


    jump_result_2 <- jump_wr_exh(X_2 / population_test_2$pi_i, population_test_2$pi_i, eps)

    it("fixes the jumping coordinate to an integer" |> add_context(k),
      expect_gt(
        sum(abs(jump_result_2 - round(jump_result_2)) < eps),
        sum(abs(population_test_2$pi_i - round(population_test_2$pi_i)) < eps)
      )
    )
    it("gives positive jumping probabilities" |> add_context(k),
      expect_true(all(jump_result_2 > -eps))
    )
    it("verifies the balancing constraints" |> add_context(k),
      expect_equal(colSums(X_2 * jump_result_2 / population_test_2$pi_i), colSums(X_2))
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

describe("flight_wr_exh()", {
  nb_aux_var_1 <- dim(X_1)[2]
  nb_aux_var_2 <- dim(X_2)[2]
  set.seed(1234)
  for (k in seq_len(50)) {
    s_1 <- flight_wr_exh(X_1, population_test_1$pi_i, eps)
    it("has at most p non-integer coordinates, where p is the number of auxiliary variables" |> add_context(k), { # nolint
      expect_lte(sum(abs(round(s_1) - s_1) > eps), nb_aux_var_1)
    })
    it("verifies the balancing constraints" |> add_context(k), {
      expect_equal(colSums(X_1 * s_1 / population_test_1$pi_i), colSums(X_1))
    })

    s_2 <- flight_wr_exh(X_2, population_test_2$pi_i, eps)
    it("has at most p non-integer coordinates, where p is the number of auxiliary variables (2)" |> add_context(k), { # nolint
      expect_lte(sum(abs(round(s_2) - s_2) > eps), nb_aux_var_2)
    })
    it("verifies the balancing constraints (2)" |> add_context(k), {
      expect_equal(colSums(X_2 * s_2 / population_test_2$pi_i), colSums(X_2))
    })
  }
})

describe("sample_cube_wr_exh()", {
  nb_aux_var_1 <- dim(X_1)[2]
  nb_aux_var_2 <- dim(X_2)[2]
  set.seed(1234)
  for (k in seq_len(50)) {
    s_1 <- sample_cube_wr_exh(X_1, population_test_1$pi_i, eps)
    it("has at most p non-integer coordinates, where p is the number of auxiliary variables" |>
        add_context(k),
      expect_lte(sum(abs(round(s_1) - s_1) > eps), nb_aux_var_1)
    )
    it("verifies the balancing constraints" |>
        add_context(k),
      expect_equal(colSums(X_1 * s_1 / population_test_1$pi_i), colSums(X_1))
    )

    s_2 <- sample_cube_wr_exh(X_2, population_test_2$pi_i, eps)
    it("has at most p non-integer coordinates, where p is the number of auxiliary variables (2)" |> add_context(k), # nolint
      expect_lte(sum(abs(round(s_2) - s_2) > eps), nb_aux_var_2)
    )
    it("verifies the balancing constraints (2)" |> add_context(k),
      expect_equal(colSums(X_2 * s_2 / population_test_2$pi_i), colSums(X_2))
    )
  }
})

describe("jump_wr_ent()", {
  set.seed(1234)
  for (k in seq_len(10)) {
    jump_result_1 <- jump_wr_ent(X_1 / population_test_1$pi_i, population_test_1$pi_i, eps)
    n_select_1 <- jump_result_1$pi_jump[jump_result_1$k_select]

    it("fixes the jumping coordinate to an integer" |> add_context(k),
      expect_true(abs(n_select_1 - round(n_select_1)) < eps)
    )
    it("gives positive jumping probabilities" |> add_context(k),
      expect_true(all(jump_result_1$pi_jump > -eps))
    )
    it("verifies the balancing constraints" |> add_context(k),
      expect_equal(colSums(X_1 * jump_result_1$pi_jump / population_test_1$pi_i), colSums(X_1))
    )

    jump_result_2 <- jump_wr_ent(X_2 / population_test_2$pi_i, population_test_2$pi_i, eps)
    n_select_2 <- jump_result_2$pi_jump[jump_result_2$k_select]

    it("fixes the jumping coordinate to an integer (2)" |> add_context(k),
      expect_true(abs(n_select_2 - round(n_select_2)) < eps)
    )
    it("gives positive jumping probabilities (2)" |> add_context(k),
      expect_true(all(jump_result_2$pi_jump > -eps))
    )
    it("verifies the balancing constraints (2)" |> add_context(k),
      expect_equal(colSums(X_2 * jump_result_2$pi_jump / population_test_2$pi_i), colSums(X_2))
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
  s_1 <- flight_wr_ent(X_1, population_test_1$pi_i, eps)

  it("verifies the balancing constraints", {
    expect_equal(colSums(X_1 * s_1 / population_test_1$pi_i), colSums(X_1))
  })


  s_2 <- flight_wr_ent(X_2, population_test_2$pi_i, eps)

  it("verifies the balancing constraints (2)", {
    expect_equal(colSums(X_2 * s_2 / population_test_2$pi_i), colSums(X_2))
  })
})
