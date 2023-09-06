create_pseudo_population <- function(sample) {
  n_sample <- nrow(sample)
  bind_rows(map(
    seq_len(n_sample),
    \(x) sample
  )) |>
    mutate(
      x = x / pi_i,
      y = y / (n_sample * pi_i),
      pi_i = 1 / n_sample
    )
}


# Tire avec remise sans équilibrer et caler sur t_x
# Deville tille - 2005

