#' Cochran's Test for Variance Heterogeneity
#' Clause 8.3.4
test_cochran <- function(cell_stats_level, alpha_straggler = 0.05, alpha_outlier = 0.01) {
  sub <- cell_stats_level[!is.na(cell_stats_level$sd), ]
  s_vals <- sub$sd
  p <- nrow(sub)
  n <- round(mean(sub$n))

  C_stat <- max(s_vals^2) / sum(s_vals^2) # Formula in 8.3.4.2
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
    statistic = C_stat,
    suspect_lab = max_lab,
    crit_5 = crit_5,
    crit_1 = crit_1,
    status = status
  )
}

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

#' Grubbs' Test for Two Outlying Means
#' Clause 8.3.5.2
test_grubbs_double <- function(cell_stats_level, alpha_straggler = 0.05, alpha_outlier = 0.01) {
  sub <- cell_stats_level[!is.na(cell_stats_level$mean), ]
  means <- sort(sub$mean)
  p <- length(means)

  s0_sq <- sum((means - mean(means))^2)

  # Two highest
  m_high <- means[1:(p - 2)]
  s_high_sq <- sum((m_high - mean(m_high))^2)
  G_high <- s_high_sq / s0_sq

  # Two lowest
  m_low <- means[3:p]
  s_low_sq <- sum((m_low - mean(m_low))^2)
  G_low <- s_low_sq / s0_sq

  G_stat <- min(G_high, G_low) # Lower value is suspicious
  type <- if (G_high < G_low) "Two Highest" else "Two Lowest"

  # One-tailed alpha = alpha/2 for two-sided Grubbs double test
  crit_5 <- grubbs_double_critical_value(p, alpha_straggler / 2)
  crit_1 <- grubbs_double_critical_value(p, alpha_outlier / 2)

  status <- "OK"
  if (G_stat < crit_1) {
    status <- "Outlier (**)"
  } else if (G_stat < crit_5) {
    status <- "Straggler (*)"
  }

  list(
    statistic = G_stat,
    type = type,
    crit_5 = crit_5,
    crit_1 = crit_1,
    status = status
  )
}
