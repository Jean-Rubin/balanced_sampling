add_context <- function(text, ...) {
  kwargs_quo <- rlang::enquos(...)
  kwargs <- list(...)
  no_names <- setdiff(seq_along(kwargs), which(names(kwargs) == ""))
  names(kwargs)[no_names] <- purrr::map_chr(kwargs_quo[no_names], rlang::quo_name)
  kwargs[] <- purrr::map(
    kwargs,
    \(x) {
      capture.output(lobstr::tree(x)) |>
        paste0(collapse = "\n")
    }
  )

  paste0(
    text, ":\n",
    purrr::imap_chr(
      kwargs,
      \(arg, name) paste0("- ", name, " = ", arg)
    ) |> paste0(collapse = "\n")
  )
}
