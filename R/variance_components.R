#' precision
#'
#' @description
#' Calculate precision variances for each level (8.6.11)
#'
#' @param clean_data Data set with original observations for each level and lab in long format, after cleaning of outliers.
#'
#' @export
precision <- function(clean_data) {
  levels <- unique(clean_data$level)
  precision_tbl <- do.call(rbind, lapply(levels, function(lvl) {
    sub_stats <- clean_data[clean_data$level == lvl, ]
    res       <- calc_variances(sub_stats) # helper function below
    res$level <- lvl
    res
  }))
}


#' calc_variances
#'
#' @description
#' Calculates the general mean (m_hat), repeatability standard deviation (s_r), between-laboratory standard deviation (s_L), and reproducibility standard deviation (s_R) using Restricted Maximum Likelihood (REML) estimation with random effects (Annex B.2).
#'
#' @param cell_data A data.frame containing the raw individual observation test results for a single level, with columns 'lab' and 'value'.
#'
#' @importFrom nlme lme
#'
#' @return A data.frame containing p, m_hat, s_r, s_L, and s_R.
#'
#' @export
calc_variances <- function(cell_data) {

  df_raw <- cell_data[!is.na(cell_data$value), ]
  p      <- length(unique(df_raw$lab))

  # fit random intercept model y ~ 1 + (1 | lab) via REML
  fit <- nlme::lme(
    fixed  = value ~ 1,
    random = ~ 1 | lab,
    method = "REML",
    data   = df_raw
  )

  # Extract variance components (s_L and s_r)
  # VarCorr(fit) returns standard deviations as character matrix
  vc      <- nlme::VarCorr(fit)
  s_L_est <- as.numeric(vc["(Intercept)", "StdDev"])
  s_r_est <- as.numeric(vc["Residual", "StdDev"])
  s_R_est <- sqrt(s_L_est^2 + s_r_est^2)

  # Extract fixed effect intercept (m_hat)
  m_hat_est <- as.numeric(nlme::fixef(fit)[1])

  res <- data.frame(
    p     = p,
    m_hat = m_hat_est,
    s_r   = s_r_est,
    s_L   = s_L_est,
    s_R   = s_R_est
  )
  return(res)
}
