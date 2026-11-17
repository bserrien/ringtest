#' grubbs_critical_value_single
#'
#' @description
#' Grubbs' Single Outlier Test Two-Tailed Critical Value (ISO Formula D.2)
#'
#' @param p number of values in the data set tested
#' @param alpha type I error (default 0.05)
#'
#' @export
grubbs_critical_value_single <- function(p, alpha = 0.05) {
  t_val <- stats::qt(1 - alpha / (2 * p), df = p - 2)
  ## Note: the ISO 5725-2:2025(en) edition contains a typo in the formula for
  ##       th critical value, typo in the parentheses in the denominator
  ## wrong: ((p - 1) * t_val) / sqrt(p * (p - 2) + t_val^2), below is correct:
  ((p - 1) * t_val) / sqrt(p * (p - 2 + t_val^2))
}
