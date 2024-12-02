# Simple analysis of bootstrap variance estimation ----------------------------
library(targets)
library(dplyr)

tar_load(population)
n_pop <- nrow(population)
n_sample <- sum(population$pi_i)

n_iter_true <- 1000L
x_names <- c("x1", "pi_i_aux")

# Sampler
sampler_base <- sampler_gen_base(x_names)
sampler_flight_exh <- sampler_gen_flight_wr_exh(x_names)
sampler_flight_ent <- sampler_gen_flight_wr_ent(x_names)

# Base sample
sample_base <- sampler_base(population)
sample_base

# True Variance from MC estimation
compute_v_trues(population, list(base = sampler_base), n_iter_true, .progress_by_sample_fn = TRUE)
# V = 9621
mc_estimate_total(population, sampler_base, n_iter_true)

## Duplication bootstrap ------------------------------------------------------
sample_wr_copy <- sample_base |>
  tidyr::uncount(n_pop / n_sample)
compute_v_trues(sample_wr_copy, list(wr_copy = sampler_base), n_iter_true, .progress_by_sample_fn = TRUE)
# V = 10671


## Exhaustion bootstrap -------------------------------------------------------
### Inflating sample_base inclusion probabilities -----------------------------
sample_base |>
  mutate(
    pi_i = pi_i * n_pop / n_sample,
    pi_i_aux = pi_i
  ) |>
  compute_v_trues(list(wr_exh = sampler_flight_exh), n_iter_true, .progress_by_sample_fn = TRUE)

# V = 52.9
# This does not work
# pi_i = 1 -> sampled with certainty

# Repeat 10 times ?
purrr::map_dbl(
  seq_len(n_iter_true),
  \(j) purrr::map_dbl(
    seq_len(n_pop / n_sample),
    \(i) {
      sample_base |>
        mutate(
          pi_i = pi_i,
          pi_i_aux = pi_i
        ) |>
        sampler_flight_exh() |>
        estimate_total()
    }
  ) |> sum(),
  .progress = TRUE
) |> var()

# V = 10907.9
# Also not very fast


## Max-entropy ----------------------------------------------------------------
### Inflating sample_base inclusion probabilities -----------------------------
sample_base |>
  mutate(
    pi_i = pi_i * n_pop / n_sample,
    pi_i_aux = pi_i
  ) |>
  compute_v_trues(
    list(wr_exh = sampler_flight_ent),
    n_iter_true, .progress_by_sample_fn = TRUE
  )
# V = ?
# About 500x times slower than the above