library(dplyr)
library(purrr)
library(ggplot2)
library(targets)
set.seed(123)

tar_source()
tar_load(population)

noises <- c(1, 2, 5, 10, 20)
x_names <- c("pi_i_aux", "x1", "x2", "x3")
n_iter_true <- 1000L
sample_fn_list_1 <- list(
  srswor = sampler_gen_srswor(),
  # base = sampler_gen_base(x_names),
  # base_flight = sampler_gen_flight_base(x_names)
  # wr_flight_exh = sampler_gen_flight_wr_exh(x_names),
  # wr_copy = sampler_gen_wr_copy(x_names),
  # wr_flight_ent = sampler_gen_flight_wr_ent(x_names),
  # wr_flight_copy = sampler_gen_flight_wr_copy(x_names)
)
v_approx_fn_list <- get_v_approx_fn_list()
fixed_sample_fn_list <- sample_fn_list_1

population_noise <- population |>
  mutate(
    y1 = y - eps + sqrt(x3) * rnorm(n(), mean = 0, sd = noises[1]),
    y2 = y - eps + sqrt(x3) * rnorm(n(), mean = 0, sd = noises[2]),
    y3 = y - eps + sqrt(x3) * rnorm(n(), mean = 0, sd = noises[3]),
    y4 = y - eps + sqrt(x3) * rnorm(n(), mean = 0, sd = noises[4]),
    y5 = y - eps + sqrt(x3) * rnorm(n(), mean = 0, sd = noises[5])
  )

population_noise <- population_noise |>
  # set_inclusion_proba(pi_gen_unif(100))
  set_inclusion_proba(pi_gen_prop_size(100, sqrt(population_noise$x3)))
  # set_inclusion_proba(pi_gen_prop_size(100, population_noise$x3))

ggplot(population_noise, aes(x = pi_i)) +
  geom_histogram(bins = 25) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Inclusion probability", y = "Number of observations")

ggsave("output/inclusion_probabilities.pdf")

# R squared 
r_squared <- purrr::map_dbl(1:5, \(i) {
  lin_reg <- lm(reformulate(c("0", x_names), paste0("y", i)), data = population_noise)
  y_i <- population_noise[[paste0("y", i)]]

  1 - mean((y_i - predict(lin_reg))^2) / mean((y_i - mean(y_i))^2)
})

# Variance approximations
diag(v_approx_fn_list$v_deville(population_noise, x_names, y_names = paste0("y", 1:5)))
diag(v_approx_fn_list$v_multi(population_noise, x_names, y_names = paste0("y", 1:5)))
diag(v_approx_fn_list$v_deville_dup(
  population_noise, x_names, y_names = paste0("y", 1:5), n_duplication = 100
))

# Cube base methods
set.seed(8945)
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

# Cube based method
set.seed(98761)
tictoc::tic()
n_iter <- 1000L
var_y_names <- paste0("y", 1:5)
pb <- txtProgressBar(min = 0, max = n_iter, initial = 0, style = 3)
y_hats_base_flight <- vector(mode = "list", n_iter)
y_hats_base <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  x <- as.matrix(population_noise[x_names])
  sample_noise <- population_noise |>
    mutate(
      s_i = sampling::fastflightcube(.env$x, pi_i, comment = FALSE),
      s_i_land = sampling::landingcube(.env$x, s_i, pi_i, comment = FALSE)
    )

  y_hats_base_flight[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()

  y_hats_base[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i_land / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()

  setTxtProgressBar(pb, i)
}
close(pb)
tictoc::toc()
readr::write_csv(bind_rows(y_hats_base_flight), "output/y_hats_base_flight.csv")
readr::write_csv(bind_rows(y_hats_base), "output/y_hats_base.csv")
bind_rows(y_hats_base_flight) |> summarize(across(everything(), var))
bind_rows(y_hats_base) |> summarize(across(everything(), var))


var_y_names <- paste0("y", 1:5)
pb <- txtProgressBar(min = 0, max = n_iter, initial = 0, style = 3)

# wr_flight_exh treated alone
set.seed(123)
sample_flight_wr_exh <- sampler_gen_flight_wr_exh(x_names)

y_hats_exh <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  sample_noise <- population_noise |>
    sample_flight_wr_exh()

  y_hats_exh[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()
  setTxtProgressBar(pb, i)
}
close(pb)
readr::write_csv(bind_rows(y_hats_exh), "output/y_hats_exh.csv")
bind_rows(y_hats_exh) |> summarize(across(everything(), var))

# wr_copy treated alone
set.seed(123)
sample_flight_wr_copy <- sampler_gen_flight_wr_copy(x_names, var_y_names)

tictoc::tic()
y_hats_copy <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  sample_noise <- population_noise |>
    sample_flight_wr_copy()

  y_hats_copy[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()
  setTxtProgressBar(pb, i)
}
close(pb)
readr::write_csv(bind_rows(y_hats_copy), "output/y_hats_copy.csv")
tictoc::toc()

# wr_ent treated alone
set.seed(123)
sample_flight_wr_ent <- sampler_gen_flight_wr_ent(x_names)

tictoc::tic()
y_hats_wr_ent <- vector(mode = "list", n_iter)
for (i in seq_len(n_iter)) {
  sample_noise <- population_noise |>
    sample_flight_wr_ent()

  y_hats_wr_ent[[i]] <- sample_noise |>
    mutate(
      across(all_of(var_y_names), \(y) y * s_i / pi_i)
    ) |>
    select(all_of(var_y_names)) |>
    as.matrix() |>
    colSums()
  setTxtProgressBar(pb, i)
}
close(pb)
readr::write_csv(bind_rows(y_hats_wr_ent), "output/y_hats_ent.csv")
tictoc::toc()

purrr::map(
  paste0("output/", c("y_hats_base_flight", "y_hats_base", "y_hats_exh", "y_hats_ent", "y_hats_copy"), ".csv"),
  \(path) {
    readr::read_csv(path) |>
      tidyr::pivot_longer(everything()) |>
      mutate(method = .env$path)
  }
) |>
  purrr::list_rbind() |>
  summarize(
    mean_y = mean(value),
    var_y = var(value),
    .by = c(name, method)
  ) |> readr::write_csv("output/results.csv")

table_result <- readr::read_csv("output/results.csv") |>
  select(name, method, var_y) |>
  mutate(method = stringr::str_extract(method, "output/y_hats_(.*).csv", group = 1)) |>
  tidyr::pivot_wider(names_from = "method", values_from = var_y)

cbind(
  bind_rows(res_list) |> select(noise, r_squared, srswor, v_multi, v_deville),
  table_result |> select(-name)
) |>
  select(noise, r_squared, srswor, base, base_flight, exh, ent, copy, v_multi, v_deville) |>
  write_table_noise_tex(
    length(fixed_sample_fn_list),
    length(v_approx_fn_list),
    "output/table/noise_complete.tex"
  )
