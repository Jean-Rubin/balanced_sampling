# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  # packages that your targets need to run
  packages = c(
    "tibble",
    "tidyr",
    "dplyr",
    "purrr",
    "sampling"
  ),
  format = "qs" # default storage format
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

tar_source()

# Repeat construction of targets
n_multi <- 10L # Number of repetitions
multi_step <- tar_map(
  unlist = FALSE,
  values = list(iter = seq_len(n_multi)),
  # Samples from cube method
  tar_target(sample_base,
    sampler_gen_base(x_names)(population)
  ),
  # Samples from with-replacement exhaustion
  tar_target(sample_wr,
    sampler_gen_flight_wr_exh(x_names)(population)
  ),
  # Pseudo population from basic cube sample
  tar_target(pseudo_population,
    create_pseudo_population(sample_base, x_names)
  ),
  # Bootstrap estimations from pseudo population
  tar_target(y_boots,
    mc_estimate_total(pseudo_population, sampler_gen_base(x_names), n_iter_boot)
  ),
  # Estimated variance from pseudo population bootstrap
  tar_target(v_hat,
    var(y_boots)
  )
)

# Table of bootstrap estimators of total:
# `id`: (non-unique) name of the form `y_boots_{i}`.
# `y_boots`: value of the bootstrap estimator.
multi_step_y_boots <- tar_combine(
  y_bootss,
  multi_step[["y_boots"]],
  command = bind_cols(!!!.x) |>
    pivot_longer(everything(), names_to = "id", values_to = "y_boots")
)

# Named vector of bootstrap variance estimators.
multi_step_v_hat <- tar_combine(
  v_hats,
  multi_step[["v_hat"]]
)

list(
  ## Parameters ----------------------------------------------------------------
  tar_target(n_tot, 1000L), # Population size
  tar_target(n_sample, 100L), # Sample size
  tar_target(n_iter_true, 500L), # Number of iteration to compute true variance
  tar_target(n_iter_boot, 100L), # Number of bootstrap iteration
  tar_target(x_names, c("x1", "pi_i_aux")), # Names of auxiliary variables
  ## Main baseline -------------------------------------------------------------
  tar_target(population, # Population
    create_population(n_tot) |>
      set_inclusion_proba(pi_gen_unif(n_sample = n_sample))
  ),
  tar_target(y_true, # Total of `y` on the population
    compute_total(population)
  ),
  tar_target(y_hats, # Estimators of total
    mc_estimate_total(population, sampler_gen_base(x_names), n_iter_true)
  ),
  tar_target(v_true, # True variance of the estimator of total
    var(y_hats)
  ),
  tar_target(y_hats_wr, # Estimators of total from with-replacement exhaustion
    mc_estimate_total(population, sampler_gen_flight_base(x_names), n_iter_true)
  ),
  tar_target(v_true_wr, # Variance associated to with-replacement exhaustion
    var(y_hats_wr)
  ),
  tar_target(v_approx_multinomial, # Multinomial approximation of variance
    compute_v_approx_multinomial(population, x_names)
  ),
  ## Single step ---------------------------------------------------------------
  tar_target(step_sample, # Basic cube sample
    sampler_gen_base(x_names)(population)
  ),
  # Basic SRSWOR sample (has an indicator of selection)
  tar_target(step_sample_srswor,
    sampler_gen_srswor()(population)
  ),
  # With-replacement exhaustion sample
  tar_target(step_sample_wr,
    sampler_gen_flight_wr_exh(x_names)(population)
  ),
  # Pseudo-population generated from basic cube sample
  tar_target(step_pseudo_population,
    create_pseudo_population(step_sample, x_names)
  ),
  # Vector of bootstrap estimators of total generated from pseudo-population
  tar_target(step_y_boots,
    mc_estimate_total(step_pseudo_population, sampler_gen_base(x_names), n_iter_boot)
  ),
  # Pseudo-population bootstrap variance estimator
  tar_target(step_v_hat,
    var(step_y_boots)
  ),
  ## Multi step ----------------------------------------------------------------
  multi_step,
  multi_step_y_boots,
  multi_step_v_hat,
  ## Analysis ------------------------------------------------------------------
  tar_target(noise_df,
    compute_noise_df(
      population,
      noises = c(1, 2, 5, 10, 20),
      x_names = c("pi_i_aux", "x1", "x2", "x3"),
      n_iter_true = 1000L,
      sample_fn_list = get_sample_fn_list,
      v_approx_fn_list = get_v_approx_fn_list()
    )
  ),
  ### Output -------------------------------------------------------------------
  tar_target(table_noise,
    write_table_noise_tex(
      noise_df,
      get_sample_fn_list(c("pi_i_aux", "x1", "x2", "x3")),
      get_v_approx_fn_list()
    ),
    format = "file"
  )
)
