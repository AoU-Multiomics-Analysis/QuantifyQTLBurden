# Utility math helpers for burden calculations.
safe_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

safe_max_abs <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(abs(x), na.rm = TRUE)
}

safe_ratio <- function(numerator, denominator) {
  numerator <- as.numeric(numerator)
  denominator <- as.numeric(denominator)
  out <- rep(NA_real_, length(numerator))
  valid <- is.finite(numerator) & is.finite(denominator) & denominator != 0
  out[valid] <- numerator[valid] / denominator[valid]
  out
}
