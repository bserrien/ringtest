#' calc_variances_iso
#'
#' @description
#' Calculates the general mean (m_hat), repeatability standard deviation (s_r), between-laboratory standard deviation (s_L), and reproducibility standard deviation (s_R) using the ISO procedure (chapter 8.4).
#'
#' @param cell_stats A data.frame containing the cell statistics of a single level; outcome of function `calc_cell_stats.`
#'
#' @return A data.frame containing p, m_hat, s_r, s_L, and s_R.
#'
#' @export
calc_variances_iso <- function(cell_stats) {
  p       <- nrow(cell_stats)
  n_i     <- cell_stats$n
  y_bar_i <- cell_stats$mean
  s_i     <- cell_stats$sd

  # Formula (3) / (4): general mean estimate: m_hat
  m_hat <- sum(n_i * y_bar_i) / sum(n_i)

  # Formula (5): repeatability variance s2_r
  s2_r <- sum((n_i - 1) * (s_i^2)) / sum(n_i - 1)

  # Formula (7): intermediate sums for s2_L
  s2_d  <- (1 / (p - 1)) * sum(n_i * (y_bar_i - m_hat)^2)
  n_bar <- (1 / (p - 1)) * (sum(n_i) - (sum(n_i^2) / sum(n_i)))

  # Formula (6): between-laboratory variance s2_L
  s2_L <- (s2_d - s2_r) / n_bar
  if (s2_L < 0) s2_L <- 0 # Clause 8.4.5.4: truncated to zero

  # Formula (8): reproducibility variance s2_R
  s2_R <- s2_r + s2_L

  res <- data.frame(
    p     = p,
    m_hat = m_hat,
    s_r   = sqrt(s2_r),
    s_L   = sqrt(s2_L),
    s_R   = sqrt(s2_R)
  )
  return(res)
}
