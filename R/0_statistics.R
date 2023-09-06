compute_total <- function(population) {
  population |>
    summarize(y_true = sum(y)) |>
    pull(y_true)
}
