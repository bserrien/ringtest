#' test_outliers
#'
#' @description
#' Detection of stragglers and outliers with Cochran's and Grubbs' tests.
#'
#' @param stats output from function calc_cell_stats
#' @param max_iter maximal number of outlier-detection-and-cleaning iterations
#'
#' @export
test_outliers <- function(stats, max_iter = 7) {

  # outlier screen per level (8.6.5)
  levels         <- unique(stats$level)
  outlier_report <- list()

  clean_stats <- stats
  for (lvl in levels) {
    sub_stats <- clean_stats[clean_stats$level == lvl, ]
    cochran   <- test_cochran(sub_stats)
    grubbs_1  <- test_grubbs_single(sub_stats)
    grubbs_2  <- test_grubbs_double(sub_stats)
    outlier_report[[as.character(lvl)]] <- list(
      cochran       = cochran,
      grubbs_single = grubbs_1,
      grubbs_double = grubbs_2
    )

    # Discard statistical outliers (**) as per Clause 8.6.9
    if (cochran$status == "Outlier (**)") {
      clean_stats$sd[clean_stats$level == lvl &
                       clean_stats$lab == cochran$suspect_lab] <- NA
    }
    if (grubbs_1$status == "Outlier (**)") {
      clean_stats$mean[clean_stats$level == lvl &
                         clean_stats$lab == grubbs_1$suspect_lab] <- NA
      clean_stats$sd[clean_stats$level == lvl &
                       clean_stats$lab == grubbs_1$suspect_lab] <- NA
    }
  }
}
