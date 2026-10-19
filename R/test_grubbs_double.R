#' test_grubbs_double
#'
#' @description
#' Grubbs' Test for Two Outlying Means (Clause 8.3.5.2)
#'
#' @param cell_stats_level output from function `calc_cell_stats` (single level)
#' @param alpha_straggler type I error for straggler detection (default .05)
#' @param alpha_outlier type I error for outlier detection (default .01)
#'
#' @export
test_grubbs_double <- function(cell_stats_level,
                               alpha_straggler = 0.05,
                               alpha_outlier = 0.01) {

  sub   <- cell_stats_level[!is.na(cell_stats_level$mean), ]
  means <- sort(sub$mean, decreasing = FALSE) # sort means in ascending order
  p     <- length(means) # number of labs
  s0_sq <- sum((means - mean(means))^2)

  # Two highest
  m_high    <- means[1:(p - 2)]
  s_high_sq <- sum((m_high - mean(m_high))^2)
  G_high    <- s_high_sq / s0_sq

  # Two lowest
  m_low    <- means[3:p]
  s_low_sq <- sum((m_low - mean(m_low))^2)
  G_low    <- s_low_sq / s0_sq

  # lower value is suspicious, minimum ratio represents maximum drop in variance
  G_stat <- min(G_high, G_low)

  # FIX: Order original sub$mean vector (unsorted) to get correct lab names
  suspect_labs_double <- if (G_high < G_low) {
    sub$lab[order(sub$mean, decreasing = TRUE)[1:2]]
  } else {
    sub$lab[order(sub$mean, decreasing = FALSE)[1:2]]
  }

  # One-tailed alpha = alpha/2 for two-sided Grubbs double test
  crit_5 <- grubbs_critical_value_double(p, alpha_straggler / 2)
  crit_1 <- grubbs_critical_value_double(p, alpha_outlier / 2)

  status <- "OK"
  if (G_stat < crit_1) {
    status <- "Outlier (**)"
  } else if (G_stat < crit_5) {
    status <- "Straggler (*)"
  }

  list(
    statistic    = G_stat,
    suspect_labs = suspect_labs_double,
    crit_5       = crit_5,
    crit_1       = crit_1,
    status       = status
  )
}
