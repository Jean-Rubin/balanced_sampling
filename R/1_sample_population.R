sample_population <- function(population) {
  population |>
    mutate(
      i = sampling::samplecube(x, pi_i, comment = FALSE)
    ) |>
    filter(i == 1) |>
    select(-i)
}

