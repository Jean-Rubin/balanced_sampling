library(dplyr)
library(purrr)
library(targets)
set.seed(42)

tar_source()
tar_load(population)

noises <- c(1, 2, 5, 10, 20)
x_names <- c("pi_i_aux", "x1", "x2", "x3")
n_iter_true <- 10000L
sample_fn_list_1 <- list(
  srswor = sampler_gen_srswor(),
  base = sampler_gen_base(x_names),
  base_flight = sampler_gen_flight_base(x_names),
  wr_flight_exh = sampler_gen_flight_wr_exh(x_names),
  # wr_copy = sampler_gen_wr_copy(x_names),
  # wr_flight_ent = sampler_gen_flight_wr_ent(x_names),
  wr_flight_copy = sampler_gen_flight_wr_copy(x_names)
)
v_approx_fn_list <- get_v_approx_fn_list()
fixed_sample_fn_list <- sample_fn_list_1

population_noise <- population |>
  mutate(
    y1 = y - eps + rnorm(n(), mean = 0, sd = noises[1]),
    y2 = y - eps + rnorm(n(), mean = 0, sd = noises[2]),
    y3 = y - eps + rnorm(n(), mean = 0, sd = noises[3]),
    y4 = y - eps + rnorm(n(), mean = 0, sd = noises[4]),
    y5 = y - eps + rnorm(n(), mean = 0, sd = noises[5])
  )


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

readr::write_csv(bind_rows(res_list), "results.csv")


# wr_ent treated alone
set.seed(42)
sample_flight_wr_ent <- sampler_gen_flight_wr_ent(x_names)
var_y_names <- paste0("y", 1:5)

n_iter <- 10000
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

readr::write_csv(bind_rows(y_hats), "y_hats.csv")
