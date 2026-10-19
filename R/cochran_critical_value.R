#' cochran_critical_value
#'
#' @description
#' Cochran's Test Critical Value (ISO Formula D.1)
#'
#' @param p number of labs
#' @param n number of observations per lab
#' @param alpha type I error (default 0.05)
#'
#' @export
cochran_critical_value <- function(p, n, alpha = 0.05) {
  df1   <- (p - 1) * (n - 1)
  df2   <- n - 1
  f_val <- stats::qf(alpha / p, df1, df2)
  1 / (1 + (p - 1) * f_val)
}
