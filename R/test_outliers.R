#' test_outliers
#'
#' @description
#' Detection of stragglers and outliers with Cochran's and Grubbs' tests.
#'
#' @param cell_stats output from function `calc_cell_stats`
#'
#' @export
test_outliers <- function(cell_stats,
                          alpha_straggler = 0.05,
                          alpha_outlier = 0.01) {

  # initialize output structure
  res               <- cell_stats[, c("lab", "level", "n")]
  res$cochran       <- NA_character_
  res$grubbs_single <- NA_character_
  res$grubbs_double <- NA_character_

  # outlier screen per level (Clause 8.6.5)
  levels <- unique(cell_stats$level)
  for (lvl in levels) {
    idx_lvl   <- which(cell_stats$level == lvl)
    sub_stats <- cell_stats[idx_lvl, ]

    cochran <- test_cochran(sub_stats)
    if (cochran$status != "OK") {
      target_row <- which(
        res$level == lvl & res$lab == cochran$suspect_lab
      )
      res$cochran[target_row] <- cochran$status
    }

    grubbs_single <- test_grubbs_single(sub_stats)
    if (grubbs_single$status != "OK") {
      target_row <- which(
        res$level == lvl & res$lab == grubbs_single$suspect_lab
      )
      res$grubbs_single[target_row] <- grubbs_single$status
    }

    # Rule: Do NOT apply double test if single test found an outlier
    # Clause 8.3.5.3
    if (grubbs_single$status == "OK") {
      grubbs_double <- test_grubbs_double(sub_stats)
      if (grubbs_double$status != "OK") {
        target_rows <- which(
          res$level == lvl & res$lab %in% grubbs_double$suspect_labs
        )
        res$grubbs_double[target_rows] <- grubbs_double$status
      }
    }
  }

  return(res)
}
