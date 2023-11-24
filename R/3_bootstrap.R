create_pseudo_population <- function(sample, x_names) {
  n_sample <- nrow(sample)
  bind_rows(map(
    seq_len(n_sample),
    \(x) sample
  )) |>
    mutate(
      across(all_of(x_names), \(x) x / pi_i),
      y = y / (.env$n_sample * pi_i),
      pi_i = 1 / .env$n_sample
    )
}


# Tire avec remise sans équilibrer et caler sur t_x
# Deville tille - 2005

