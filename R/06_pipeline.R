#' Master Automated ISO 5725-2 Analysis Workflow
#'
#' @param data Data frame with columns 'lab', 'level', 'value'
#' @return A list containing intermediate cell statistics, outlier analysis results, level-by-level precision parameters, and fitted functional models.
run_iso5725_analysis <- function(data) {
  # 1. Compute initial cell means and standard deviations (8.6.3)
  stats <- calc_cell_stats(data)
  stats <- calc_mandel_hk(stats)

  # 2. Outlier Screen per Level (8.6.5)
  levels <- unique(stats$level)
  outlier_report <- list()

  clean_stats <- stats
  for (lvl in levels) {
    sub_stats <- clean_stats[clean_stats$level == lvl, ]

    cochran <- test_cochran(sub_stats)
    grubbs_1 <- test_grubbs_single(sub_stats)
    grubbs_2 <- test_grubbs_double(sub_stats)

    outlier_report[[as.character(lvl)]] <- list(
      cochran = cochran,
      grubbs_single = grubbs_1,
      grubbs_double = grubbs_2
    )

    # Discard statistical outliers (**) as per Clause 8.6.9
    if (cochran$status == "Outlier (**)") {
      clean_stats$sd[clean_stats$level == lvl & clean_stats$lab == cochran$suspect_lab] <- NA
    }
    if (grubbs_1$status == "Outlier (**)") {
      clean_stats$mean[clean_stats$level == lvl & clean_stats$lab == grubbs_1$suspect_lab] <- NA
      clean_stats$sd[clean_stats$level == lvl & clean_stats$lab == grubbs_1$suspect_lab] <- NA
    }
  }

  # 3. Calculate Precision Variances for Each Level (8.6.11)
  precision_tbl <- do.call(rbind, lapply(levels, function(lvl) {
    sub_stats <- clean_stats[clean_stats$level == lvl, ]
    res <- calc_iso_variances(sub_stats)
    res$level <- lvl
    res
  }))

  # 4. Fit Functional Relationships s_r(m) and s_R(m) (8.6.13)
  rel_r <- fit_precision_relationship(precision_tbl$m_hat, precision_tbl$s_r, rel = "I")
  rel_R <- fit_precision_relationship(precision_tbl$m_hat, precision_tbl$s_R, rel = "II")

  list(
    cell_statistics = stats,
    clean_statistics = clean_stats,
    outliers = outlier_report,
    level_precision = precision_tbl,
    model_repeatability = rel_r,
    model_reproducibility = rel_R
  )
}
