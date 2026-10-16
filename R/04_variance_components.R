#' Compute ISO 5725-2 Precision Parameters for a Single Level
#' Clause 8.4
calc_iso_variances <- function(cell_stats_level, method = c("ISO", "ANOVA")) {
  method <- match.arg(method)
  sub <- cell_stats_level[!is.na(cell_stats_level$mean), ]

  p <- nrow(sub)
  n_i <- sub$n
  y_bar_i <- sub$mean
  s_i <- sub$sd

  # Formula (3): General Mean estimate m_hat
  m_hat <- sum(n_i * y_bar_i) / sum(n_i)

  if (method == "ISO") {
    # Formula (5): Repeatability Variance s_r^2
    s_r_sq <- sum((n_i - 1) * (s_i^2), na.rm = TRUE) / sum(n_i - 1)

    # Formula (7): Intermediate sums for s_L^2
    s_d_sq <- (1 / (p - 1)) * sum(n_i * (y_bar_i - m_hat)^2)
    n_bar <- (1 / (p - 1)) * (sum(n_i) - (sum(n_i^2) / sum(n_i)))

    # Formula (6): Between-Laboratory Variance s_L^2
    s_L_sq <- (s_d_sq - s_r_sq) / n_bar
    if (s_L_sq < 0) s_L_sq <- 0 # Clause 8.4.5.4: truncated to zero

  } else if (method == "ANOVA") {
    # One-Way ANOVA approach (Annex B.1)
    # Balanced check
    n <- n_i[1]
    SS_b <- n * sum((y_bar_i - m_hat)^2)
    SS_w <- sum((n_i - 1) * s_i^2)
    MS_b <- SS_b / (p - 1)
    MS_w <- SS_w / (p * (n - 1))

    s_r_sq <- MS_w # Formula (B.1)
    s_L_sq <- (MS_b - MS_w) / n # Formula (B.2)
    if (s_L_sq < 0) s_L_sq <- 0
  }

  # Formula (8): Reproducibility Variance s_R^2
  s_R_sq <- s_r_sq + s_L_sq

  data.frame(
    p = p,
    m_hat = m_hat,
    s_r = sqrt(s_r_sq),
    s_L = sqrt(s_L_sq),
    s_R = sqrt(s_R_sq)
  )
}
