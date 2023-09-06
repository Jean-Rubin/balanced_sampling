estimate_total <- function(sample) {
  sample |>
    summarize(y_hat = sum(y / pi_i)) |>
    pull(y_hat)
}
