library(targets)
library(dplyr)
tar_source()
tar_load(population)
n_iter_true <- 1000L
x_names <- c("pi_i_aux", "const", "x1", "x2", "x3")
n_samples <- c(50, 100, 200, 500, 750)

sample_fn_list <- list(
  base = sampler_gen_base(x_names),
  # wr_copy = sampler_gen_wr_copy(x_names), # very costly
  srswor = sampler_gen_srswor()
)

v_approx_fn_list <- list(
  v_multi = compute_v_approx_multinomial,
  v_deville = compute_v_approx_deville_tille
)

sample_flight_fn_list <- list(
  wr_flight_exh = sampler_gen_flight_wr_exh(x_names)
)

# Check Baseline simulation for different sample size --------------------------
design_df <- tibble(
  n_sample = n_samples
) |>
  mutate(
    stats = purrr::map(
      n_sample,
      \(n_s) {
        pi_gen <- pi_gen_unif(n_s)
        population <- set_inclusion_proba(population, pi_gen)

        c(
          compute_v_trues(population, sample_fn_list, n_iter_true),
          compute_v_approxs(population, v_approx_fn_list, x_names)
        )
      },
      .progress = TRUE
    )
  ) |>
  tidyr::unnest_wider(stats)

# n_iter_true : 1000
# n_sample   base  srswor v_multi v_deville
# <dbl>  <dbl>   <dbl>   <dbl>     <dbl>
# 1       50 26837. 337602.  20540.    19552.
# 2      100 10866. 157641.  10270.     9261.
# 3      200  4318.  66037.   5135.     4116.
# 4      500  1099.  18760.   2054.     1029.
# 5      750   374.   6248.   1369.      343.




# Flight phase analysis --------------------------------------------------------

flight_design_df <- tibble(
  n_sample = c(50, 100, 200, 300, 400, 500, 750)
) |>
  mutate(
    stats = purrr::map(
      n_sample,
      \(n_s) {
        pi_gen <- pi_gen_unif(n_s)
        population <- set_inclusion_proba(population, pi_gen)

        c(
          compute_v_trues(population, sample_flight_fn_list, n_iter_true),
          compute_v_approxs(population, v_approx_fn_list, x_names)
        )
      },
      .progress = TRUE
    )
  ) |>
  tidyr::unnest_wider(stats)

# n_iter_true : 1000
# n_sample wr_fuzzy v_multi v_deville
# <dbl>    <dbl>   <dbl>     <dbl>
# 1       50   21932.  20540.    19552.
# 2      100   11257.  10270.     9261.
# 3      200    5583.   5135.     4116.
# 4      300    3265.   3423.     2401.
# 5      400    2840.   2567.     1544.
# 6      500    1910.   2054.     1029.
# 7      750    1304.   1369.      343.


# n_iter_true : 1000
#   n_sample wr_flight_exh v_multi v_deville
#      <dbl>         <dbl>   <dbl>     <dbl>
# 1       50        21130.  20740.    19802.
# 2      100        10774.  10370.     9380.
# 3      200         5284.   5185.     4169.
# 4      300         3355.   3457.     2432.
# 5      400         2451.   2592.     1563.
# 6      500         1883.   2074.     1042.
# 7      750         1065.   1383.      347.