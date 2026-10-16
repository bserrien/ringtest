#' Exact Critical Values and Indicators for ISO 5725-2 (Annex D)

#' Cochran's Test Critical Value
#' Formula D.1
cochran_critical_value <- function(p, n, alpha = 0.05) {
  df1 <- (p - 1) * (n - 1)
  df2 <- n - 1
  f_val <- stats::qf(1 - (alpha / p), df1, df2)
  1 / (1 + (p - 1) * f_val)
}

#' Grubbs' Single Outlier Test Two-Tailed Critical Value
#' Formula D.2
grubbs_single_critical_value <- function(p, alpha = 0.05) {
  t_val <- stats::qt(1 - alpha / (2 * p), df = p - 2)
  ((p - 1) * t_val) / sqrt(p * (p - 2 + t_val^2))
}

#' Grubbs' Double Outlier Test One-Tailed Critical Value
#' Formulae D.3 & D.4
grubbs_double_critical_value <- function(p, alpha = 0.05) {
  # Coefficients table from Table D.1
  coeff_tbl <- data.frame(
    a = c(0.001, 0.005, 0.010, 0.025, 0.050, 0.100),
    g0 = c(-4.2493, -3.6613, -3.3101, -2.8580, -2.5075, -2.1615),
    g1 = c(1.0012, 0.9558, 0.9250, 0.8833, 0.8501, 0.8169),
    g2 = c(0.0443, 0.0388, 0.0362, 0.0322, 0.0289, 0.0251)
  )
  idx <- which.min(abs(coeff_tbl$a - alpha))
  g0 <- coeff_tbl$g0[idx]
  g1 <- coeff_tbl$g1[idx]
  g2 <- coeff_tbl$g2[idx]

  f_n <- g0 + g1 * p + g2 * p^2
  f_val <- stats::qf((1 - alpha)^(1 / f_n), df1 = 2, df2 = p - 3)
  1 / (1 + (2 / (p - 3)) * f_val)
}

#' Mandel's h Indicator Line
#' Formula D.5
mandel_h_indicator <- function(p, alpha = 0.05) {
  t_val <- stats::qt(1 - alpha / 2, df = p - 2)
  ((p - 1) * t_val) / sqrt(p * (p - 2 + t_val^2))
}

#' Mandel's k Indicator Line
#' Formula D.6
mandel_k_indicator <- function(p, n, alpha = 0.05) {
  df1 <- (p - 1) * (n - 1)
  df2 <- n - 1
  f_val <- stats::qf(alpha, df1, df2)
  sqrt(p / (1 + (p - 1) * f_val))
}
