#' grubbs_critical_value_double
#'
#' @description
#' Grubbs' Double Outlier Test One-Tailed Critical Value (ISO Formulae D.3-.4)
#'
#' @param p number of labs
#' @param alpha type I error should be one of 0.001, 0.005, 0.01, 0.025, 0.05 (default), 0.1
#'
#' @export
grubbs_critical_value_double <- function(p, alpha = 0.05) {
  stopifnot(
    "alpha should be one of 0.001, 0.005, 0.01, 0.025, 0.05 (default), 0.1" =
    alpha %in% c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1)
  )

  # Coefficients table from Table D.1
  coeff_tbl <- data.frame(
    a  = c(0.001, 0.005, 0.010, 0.025, 0.050, 0.100),
    g0 = c(-4.2493, -3.6613, -3.3101, -2.8580, -2.5075, -2.1615),
    g1 = c(1.0012, 0.9558, 0.9250, 0.8833, 0.8501, 0.8169),
    g2 = c(0.0443, 0.0388, 0.0362, 0.0322, 0.0289, 0.0251)
  )
  idx   <- which.min(abs(coeff_tbl$a - alpha))
  g0    <- coeff_tbl$g0[idx]
  g1    <- coeff_tbl$g1[idx]
  g2    <- coeff_tbl$g2[idx]
  f_n   <- g0 + g1 * p + g2 * p^2
  f_val <- stats::qf((1 - alpha)^(1 / f_n), df1 = 2, df2 = p - 3)
  1 / (1 + (2 / (p - 3)) * f_val)
}
