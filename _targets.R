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

n_multi <- 10L

multi_step <- tar_map(
  unlist = FALSE,
  values = list(iter = seq_len(n_multi)),
  tar_target(sample_base,
    sampler_gen(x_names)(population)
  ),
  tar_target(sample_wr,
    sampler_gen_fuzzy_wr(x_names)(population)
  ),
  tar_target(pseudo_population,
    create_pseudo_population(sample_base, x_names)
  ),
  tar_target(y_boots,
    mcmc_population(pseudo_population, sampler_gen(x_names), n_iter_boot)
  ),
  tar_target(v_hat,
    var(y_boots)
  )
)

multi_step_y_boots <- tar_combine(
  y_bootss,
  multi_step[["y_boots"]],
  command = bind_cols(!!!.x) |>
    pivot_longer(everything(), names_to = "id", values_to = "y_boots")
)

multi_step_v_hat <- tar_combine(
  v_hats,
  multi_step[["v_hat"]]
)

list(
  ## Parameters ----------------------------------------------------------------
  tar_target(n_tot, 1000L),
  tar_target(n_sample, 100L),
  tar_target(n_iter_true, 500L),
  tar_target(n_iter_boot, 100L),
  tar_target(x_names, c("x1", "pi_i_aux")),
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
  tar_target(y_hats_wr,
    mcmc_population(population, sampler_gen_fuzzy_wr(x_names), n_iter_true)
  ),
  tar_target(v_true_wr,
    var(y_hats_wr)
  ),
  tar_target(v_approx_multinomial,
    compute_v_approx_multinomial(population, x_names)
  ),
  ## Single step ---------------------------------------------------------------
  tar_target(step_sample,
    sampler_gen(x_names)(population)
  ),
  tar_target(step_sample_srswor,
    sampler_gen_srswor()(population)
  ),
  tar_target(step_sample_wr,
    sampler_gen_fuzzy_wr(x_names)(population)
  ),
  tar_target(step_pseudo_population,
    create_pseudo_population(step_sample, x_names)
  ),
  tar_target(step_y_boots,
    mcmc_population(step_pseudo_population, sampler_gen(x_names), n_iter_boot)
  ),
  tar_target(step_v_hat,
    var(step_y_boots)
  ),
  ## Multi step ----------------------------------------------------------------
  multi_step,
  multi_step_y_boots,
  multi_step_v_hat
)
