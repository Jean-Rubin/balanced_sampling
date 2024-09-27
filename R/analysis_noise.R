# Computing variance for varying level of noises ------------------------------
# Variance estimation and variance approximations are done

# Also:
# Checking that flight phase variance weight more than landing phase
# if correlation is weaker

#' Get list of sample strategies
#'
#' These are mostly balanced sampling algorithms
#'
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#'
#' @return A named list with sampling functions.
#' @export
#' @examples
#' get_sample_fn_list(c("const", "pi_i_aux"))
get_sample_fn_list <- function(x_names) {
  list(
    srswor = sampler_gen_srswor(),
    base = sampler_gen(x_names),
    base_flight = sampler_gen_fuzzy_custom(x_names),
    wr_flight = sampler_gen_fuzzy_wr(x_names),
    wr_copy = sampler_gen_wr_copy(x_names),
    wr_ent_flight = sampler_gen_fuzzy_wr_ent(x_names)
  )
}

#' Get list of variance approximations
#'
#' @return A named list with variance approximations.
#' @export
#' @examples
#' get_v_approx_fn_list()
get_v_approx_fn_list <- function() {
  list(
    v_multi = compute_v_approx_multinomial,
    v_deville = compute_v_approx_deville_tille
  )
}

#' Generate a table comparing variance under different levels of noise
#'
#' @param population A data frame corresponding to the population.
#' @param noises A vector with the varying levels of noise.
#' @param x_names A vector of the names of the auxiliary variables used to
#'   balance the sample.
#' @param n_iter_true Number of iterations to compute variance by Monte Carlo
#'   estimation.
#' @param sample_fn_list A function returning a named list of sampling
#'   functions. The function takes `x_names` as an argument.
#' @param v_approx_fn_list A named list of variance approximations.
#'
#' @return A table with the computed variance for each level of noise.
#' @export
#'
#' @examples
#' population <- data.frame(y = 1:3, x1 = 2:4, pi_i = c(0.2, 0.3, 0.5))
#' compute_noise_df(
#'   population,
#'   noises = c(0.1, 0.5, 1, 2, 5, 10),
#'   x_names = "x1",
#'   n_iter_true = 100L,
#'   sample_fn_list = \(x_names) list(base = sampler_gen(x_names)),
#'   v_approx_fn_list = list(v_multi = compute_v_approx_multinomial)
#' )
compute_noise_df <- function(
  population, noises, x_names, n_iter_true,
  sample_fn_list, v_approx_fn_list
) {
  fixed_sample_fn_list <- sample_fn_list(x_names)

  tibble(
    noise = noises
  ) |>
    mutate(stats = purrr::map(
      noise,
      \(noise) {
        population_noise <- population |>
          mutate(
            # TODO expression is hard coded
            lin_rel = 2 * x1,
            y = lin_rel + rnorm(n(), mean = 0, sd = noise)
          )

        r_squared <- population_noise |>
          summarize(r_sq = 1 - mean((y - lin_rel)^2) / mean((y - mean(y))^2)) |>
          pull(r_sq)

        c(
          r_squared = r_squared,
          compute_v_trues(population_noise, fixed_sample_fn_list, n_iter_true),
          compute_v_approxs(population_noise, v_approx_fn_list, x_names)
        )
      },
      .progress = FALSE
    )) |>
    tidyr::unnest_wider(stats)
}

#' Write the table of noise
#'
#' The table is written in the LaTeX format.
#'
#' @param result_df Table result of the experiment for varying level of noise.
#' @param fixed_sample_fn_list List of sampling functions that where used for
#'   producing the table.
#' @param v_approx_fn_list List of approximation functions that where used for
#'   producing the table.
#'
#' @return The path where the table was saved.
#' @export
write_table_noise_tex <- function(
  result_df, fixed_sample_fn_list, v_approx_fn_list
) {
  result_df |>
    mutate(across(-noise, \(x) signif(x, digits = 3))) |>
    kableExtra::kbl(
      format = "latex",
      booktabs = TRUE,
      caption = "Influence of noise on variance",
      format.args = list(digits = 3, big.mark = ",")
    ) |>
    kableExtra::add_header_above(c(
      "", "",
      "Variance MC" = length(fixed_sample_fn_list),
      "Variance Approximation" = length(v_approx_fn_list)
    )) |>
    kableExtra::save_kable("output/table/noise.tex")

  "output/table/noise.tex"
}
