#' test_grubbs_single
#'
#' @description
#' Grubbs' test for single outlying mean (Clause 8.3.5.1)
#'
#' @param cell_stats_level output from function `calc_cell_stats` (single level)
#' @param alpha_straggler type I error for straggler detection (default .05)
#' @param alpha_outlier type I error for outlier detection (default .01)
#'
#' @export
test_grubbs_single <- function(cell_stats_level,
                               alpha_straggler = 0.05,
                               alpha_outlier = 0.01) {

  sub   <- cell_stats_level[!is.na(cell_stats_level$mean), ]
  means <- sub$mean
  p     <- length(means)    # number of labs
  x_bar <- mean(means)      # mean across labs
  s     <- stats::sd(means) # sd across labs

  G_high <- (max(means) - x_bar) / s
  G_low  <- (x_bar - min(means)) / s
  G_stat <- max(G_high, G_low)

  suspect_lab <- if (G_high > G_low) {
    sub$lab[which.max(means)]
  } else {
    sub$lab[which.min(means)]
  }

  crit_5 <- grubbs_critical_value_single(p, alpha_straggler)
  crit_1 <- grubbs_critical_value_single(p, alpha_outlier)

  status <- "OK"
  if (G_stat > crit_1) {
    status <- "Outlier (**)"
  } else if (G_stat > crit_5) {
    status <- "Straggler (*)"
  }

  list(
    statistic   = G_stat,
    suspect_lab = suspect_lab,
    crit_5      = crit_5,
    crit_1      = crit_1,
    status      = status
  )
}
