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

n_multi <- 100L

multi_step <- tar_map(
  unlist = FALSE,
  values = list(iter = seq_len(n_multi)),
  tar_target(sample,
    sampler_gen(x_names)(population)
  ),
  tar_target(sample_copy,
    sampler_gen_wr_copy(x_names)(population)
  ),
  tar_target(pseudo_population,
    create_pseudo_population(sample)
  ),
  tar_target(y_boots,
    mcmc_population(pseudo_population, sampler_gen(x_names), n_iter_boot)
  ),
  tar_target(v_hat,
    var(y_boots)
  )
)

y_boots_m <- tar_combine(
  y_bootss,
  multi_step[["y_boots"]],
  command = bind_cols(!!!.x) |>
    pivot_longer(everything(), names_to = "id", values_to = "y_boots")
)

v_hat_m <- tar_combine(
  v_hats,
  multi_step[["v_hat"]]
)

list(
  ## Parameters ----------------------------------------------------------------
  tar_target(n_tot, 1000L),
  tar_target(n_sample, 100L),
  tar_target(n_iter_true, 50L),
  tar_target(n_iter_boot, 10L),
  tar_target(x_names, c("x", "pi_i_aux")),
  ## Main baseline -------------------------------------------------------------
  tar_target(population,
    create_population(n_tot) |>
      set_inclusion_proba(pi_gen_unif(n_sample = n_sample))
  ),
  tar_target(y_true,
    compute_total(population)
  ),
  tar_target(y_hats,
    mcmc_population(population, sampler_gen(x_names), n_iter_true)
  ),
  tar_target(v_true,
    var(y_hats)
  ),
  tar_target(y_hats_wr_copy,
    mcmc_population(population, sampler_gen_wr_copy(x_names), n_iter_true)
  ),
  tar_target(v_true_wr_copy,
    var(y_hats_wr_copy)
  ),
  tar_target(v_approx_multinomial,
    compute_v_approx_multinomial(population, c("x", "pi_i_aux"))
  ),
  ## Single step ---------------------------------------------------------------
  tar_target(sample_step,
    sampler_gen(x_names)(population)
  ),
  tar_target(sample_step_srswor,
    sampler_gen_srswor()(population)
  ),
  tar_target(sample_step_wr_copy,
    sampler_gen_wr_copy(x_names)(population)
  ),
  tar_target(pseudo_population_step,
    create_pseudo_population(sample_step)
  ),
  tar_target(y_boots_step,
    mcmc_population(pseudo_population_step, sampler_gen(x_names), n_iter_boot)
  ),
  tar_target(v_hat_step,
    var(y_boots_step)
  ),
  ## Multi step ----------------------------------------------------------------
  multi_step,
  y_boots_m,
  v_hat_m
)
