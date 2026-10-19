#' test_cochran
#'
#' @description
#' Cochran's Test for Variance Heterogeneity (Clause 8.3.4)
#'
#' @param cell_stats_level output from function `calc_cell_stats` (single level)
#' @param alpha_straggler type I error for straggler detection (default .05)
#' @param alpha_outlier type I error for outlier detection (default .01)
#'
#' @export
test_cochran <- function(cell_stats_level,
                         alpha_straggler = 0.05,
                         alpha_outlier = 0.01) {

  sub    <- cell_stats_level[!is.na(cell_stats_level$sd), ]
  s_vals <- sub$sd
  p      <- nrow(sub)          # number of labs
  n      <- round(mean(sub$n)) # average number of observations per level
                               # see 8.3.4.3 for explanation, caveat

  C_stat  <- max(s_vals^2) / sum(s_vals^2) # Formula in 8.3.4.2
  max_lab <- sub$lab[which.max(s_vals^2)]

  crit_5 <- cochran_critical_value(p, n, alpha_straggler)
  crit_1 <- cochran_critical_value(p, n, alpha_outlier)

  status <- "OK"
  if (C_stat > crit_1) {
    status <- "Outlier (**)"
  } else if (C_stat > crit_5) {
    status <- "Straggler (*)"
  }

  list(
    statistic   = C_stat,
    suspect_lab = max_lab,
    crit_5      = crit_5,
    crit_1      = crit_1,
    status      = status
  )
}
