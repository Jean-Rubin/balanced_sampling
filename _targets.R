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


multi_step <- tar_map(
  unlist = FALSE,
  values = list(iter = seq_len(1000)),
  tar_target(sample,
    sample_population(population)
  ),
  tar_target(pseudo_population,
    create_pseudo_population(sample)
  ),
  tar_target(y_boots,
    mcmc_population(pseudo_population, n_iter_true)
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
  tar_target(n_iter_true, 10000L),
  ## Main baseline -------------------------------------------------------------
  tar_target(population,
    create_population(n_tot)
  ),
  tar_target(y_true,
    compute_total(population)
  ),
  tar_target(y_hats,
    mcmc_population(population, n_iter_true)
  ),
  tar_target(v_true,
    var(y_hats)
  ),
  ## Single step ---------------------------------------------------------------
  tar_target(sample_step,
    sample_population(population)
  ),
  tar_target(pseudo_population_step,
    create_pseudo_population(sample_step)
  ),
  tar_target(y_boots_step,
    mcmc_population(pseudo_population_step, n_iter_true)
  ),
  tar_target(v_hat_step,
    var(y_boots_step)
  ),
  ## Multi step ----------------------------------------------------------------
  multi_step,
  y_boots_m,
  v_hat_m
)

