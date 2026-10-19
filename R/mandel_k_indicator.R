#' mandel_k_indicator
#'
#' @description
#' Mandel's k Indicator Line (ISO Formula D.6)
#'
#' @param p number of labs
#' @param n number of observations in each lab
#' @param alpha type I error, default 0.05
#'
#' @export
mandel_k_indicator <- function(p, n, alpha = 0.05) {
  df1   <- (p - 1) * (n - 1)
  df2   <- n - 1
  f_val <- stats::qf(alpha, df1, df2)
  sqrt(p / (1 + (p - 1) * f_val))
}
