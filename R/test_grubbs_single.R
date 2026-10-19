#' Grubbs' Test for Single Outlying Mean
#' Clause 8.3.5.1
test_grubbs_single <- function(cell_stats_level, alpha_straggler = 0.05, alpha_outlier = 0.01) {
  sub <- cell_stats_level[!is.na(cell_stats_level$mean), ]
  means <- sub$mean
  p <- length(means)
  x_bar <- mean(means)
  s <- stats::sd(means)

  G_high <- (max(means) - x_bar) / s
  G_low  <- (x_bar - min(means)) / s

  G_stat <- max(G_high, G_low)
  suspect_lab <- if (G_high > G_low) sub$lab[which.max(means)] else sub$lab[which.min(means)]

  crit_5 <- grubbs_single_critical_value(p, alpha_straggler)
  crit_1 <- grubbs_single_critical_value(p, alpha_outlier)

  status <- "OK"
  if (G_stat > crit_1) {
    status <- "Outlier (**)"
  } else if (G_stat > crit_5) {
    status <- "Straggler (*)"
  }

  list(
    statistic = G_stat,
    suspect_lab = suspect_lab,
    crit_5 = crit_5,
    crit_1 = crit_1,
    status = status
  )
}
