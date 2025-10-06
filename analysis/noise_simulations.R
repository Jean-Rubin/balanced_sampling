library(dplyr)
library(purrr)
library(targets)
set.seed(123)

tar_source()
tar_load(population)

noises <- c(1, 2, 5, 10, 20)
x_names <- c("pi_i_aux", "x1", "x2", "x3")
n_iter_true <- 10000L
sample_fn_list_1 <- list(
  srswor = sampler_gen_srswor(),
  base = sampler_gen_base(x_names),
  base_flight = sampler_gen_flight_base(x_names)
  # wr_flight_exh = sampler_gen_flight_wr_exh(x_names),
  # wr_copy = sampler_gen_wr_copy(x_names),
  # wr_flight_ent = sampler_gen_flight_wr_ent(x_names),
  # wr_flight_copy = sampler_gen_flight_wr_copy(x_names)
)
v_approx_fn_list <- get_v_approx_fn_list()
fixed_sample_fn_list <- sample_fn_list_1

set.seed(123)
population_noise <- population |>
  mutate(
    y1 = y - eps + rnorm(n(), mean = 0, sd = noises[1]),
    y2 = y - eps + rnorm(n(), mean = 0, sd = noises[2]),
    y3 = y - eps + rnorm(n(), mean = 0, sd = noises[3]),
    y4 = y - eps + rnorm(n(), mean = 0, sd = noises[4]),
    y5 = y - eps + rnorm(n(), mean = 0, sd = noises[5])
  )

# R squared 
r_squared <- purrr::map_dbl(1:5, \(i) {
  lin_reg <- lm(reformulate(c("0", x_names), paste0("y", i)), data = population_noise)
  y_i <- population_noise[[paste0("y", i)]]
  
  1 - mean((y_i - predict(lin_reg))^2) / mean((y_i - mean(y_i))^2)
})

# Variance approximations
diag(v_approx_fn_list$v_deville(population_noise, x_names, y_names = paste0("y", 1:5)))
diag(v_approx_fn_list$v_multi(population_noise, x_names, y_names = paste0("y", 1:5)))


# Cube base methods
set.seed(123)
res_list <- vector(mode = "list", length = length(noises))
pb <- txtProgressBar(min = 0, max = length(noises), initial = 0, style = 3)
for (i in seq_along(noises)) {
  noise <- noises[i]
  population_noise[, "y"] <- population_noise[, paste0("y", i)]

  lin_reg <- lm(reformulate(c("0", x_names), "y"), data = population_noise)
  r_squared <- population_noise |>
    summarize(r_sq = 1 - mean((y - predict(lin_reg))^2) / mean((y - mean(y))^2)) |>
    pull(r_sq)

  res_list[[i]] <- c(
    noise = noise,
    r_squared = r_squared,
    # Should be faster to adapt and use the same simulations for multiple variables of interest
    compute_v_trues(population_noise, fixed_sample_fn_list, n_iter_true),
    compute_v_approxs(population_noise, v_approx_fn_list, x_names)
  )
  setTxtProgressBar(pb, i)
}
close(pb)
bind_rows(res_list)

write_table_noise_tex(
  bind_rows(res_list),
  length(fixed_sample_fn_list),
  length(v_approx_fn_list),
  "output/table/noise.tex"
)


n_iter <- 10000L
var_y_names <- paste0("y", 1:5)
pb <- txtProgressBar(min = 0, max = n_iter, initial = 0, style = 3)

# wr_flight_exh treated alone
set.seed(123)
sample_flight_wr_exh <- sampler_gen_flight_wr_exh(x_names)

y_hats <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  sample_noise <- population_noise |>
    sample_flight_wr_exh()

  y_hats[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()
  setTxtProgressBar(pb, i)
}
close(pb)
readr::write_csv(bind_rows(y_hats), "output/y_hats_exh.csv")
bind_rows(y_hats) |> summarize(across(everything(), var))

# wr_copy treated alone
set.seed(123)
sample_flight_wr_copy <- sampler_gen_flight_wr_copy(x_names, var_y_names)

tictoc::tic()
y_hats <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  sample_noise <- population_noise |>
    sample_flight_wr_copy()

  y_hats[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()
  setTxtProgressBar(pb, i)
}
close(pb)
readr::write_csv(bind_rows(y_hats), "output/y_hats_copy.csv")
tictoc::toc()

# wr_ent treated alone
set.seed(123)
sample_flight_wr_ent <- sampler_gen_flight_wr_ent(x_names)

tictoc::tic()
y_hats <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  sample_noise <- population_noise |>
    sample_flight_wr_ent()

  y_hats[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()
  setTxtProgressBar(pb, i)
}
close(pb)
readr::write_csv(bind_rows(y_hats), "output/y_hats_ent.csv")
tictoc::toc()

purrr::map(
  paste0("output/", c("y_hats_exh", "y_hats_ent", "y_hats_copy"), ".csv"),
  \(path) {
    readr::read_csv(path) |>
      summarize(across(everything(), var))
  }
) |>
  purrr::list_rbind()

readr::read_csv("output/results.csv")
