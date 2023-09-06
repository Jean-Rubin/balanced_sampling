compute_relative_bias <- function(v_hats, v_true) {
  (mean(v_hats) - v_true) / v_true
}

compute_relative_stability <- function(v_hats, v_true) {
  sqrt(mean((v_hats - v_true)^2)) / v_true
}

