#' outlier_removal_iso
#'
#' @description
#' Iterative process to detect and remove outliers as described in the ISO 5725-2 manual.
#'
#' @param cell_stats A data.frame containing the cell statistics; outcome of function `calc_cell_stats.`
#' @param maxiter Maximal number of iterations to perform outlier detection and removal.
#'
#' @return A data.frame containing the cell statistics of labs and levels after outlier removal. Intermediate steps are printed to the console when `verbose` = TRUE (default).
#'
#' @export
outlier_removal_iso <- function(cell_stats, maxiter = NULL, verbose = TRUE) {

  # new data.frame to perform outlier detection & cleaning on
  cell_stats_clean <- cell_stats

  # take (nrow - 1) as max number of iterations if not provided
  if (is.null(maxiter)) maxiter <- nrow(cell_stats) - 1

  for (iter in 1:maxiter) {

    # print outliers & stragglers
    outl_strag_i <- cell_stats_clean |>
      test_outliers() |>
      subset(
        !is.na(cochran) |
          !is.na(grubbs_single) |
          !is.na(grubbs_double)
      )
    if (verbose) {
      cat("iteration ", iter, ":\n", sep = "")
      print(outl_strag_1)
    }

    outl_i <- outl_strag_i |>
      subset(
        cochran == "Outlier (**)" |
          grubbs_single == "Outlier (**)" |
          grubbs_double == "Outlier (**)"
      )
    # stop if no more outliers identified
    if (nrow(outl_i) == 0) {
      break
    } else {
      # remove outliers, retain stragglers
      cell_stats_clean <- dplyr::anti_join(
        cell_stats_clean, outl_i, by = c("lab", "level")
      )
    }
  }

  return(cell_stats_clean)
}
