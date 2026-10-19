#' mandel_h_indicator
#'
#' @description
#' Mandel's h Indicator Line (ISO Formula D.5)
#'
#' @param p number of labs
#' @param alpha type I error, default 0.05
#'
#' @export
mandel_h_indicator <- function(p, alpha = 0.05) {
  t_val <- stats::qt(1 - alpha / 2, df = p - 2)
  ((p - 1) * t_val) / sqrt(p * (p - 2 + t_val^2))
}
