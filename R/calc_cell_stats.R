#' calc_cell_stats
#'
#' @description
#' Calculate cell means and within-cell standard deviations
#'
#' @param data A data.frame with columns: `lab`, `level`, `value`
#'
#' @export
calc_cell_stats <- function(data) {
  labs   <- unique(data$lab)
  levels <- unique(data$level)

  res <- expand.grid(lab = labs, level = levels, stringsAsFactors = FALSE)
  res$n    <- 0
  res$mean <- NA_real_
  res$sd   <- NA_real_

  for (i in seq_len(nrow(res))) {
    sub <- data[data$lab == res$lab[i] & data$level == res$level[i], "value"]
    sub <- sub[!is.na(sub)]
    res$n[i] <- length(sub)
    if (length(sub) > 0) {
      res$mean[i] <- mean(sub)
    }
    if (length(sub) > 1) {
      res$sd[i] <- stats::sd(sub)
    }
  }
  return(res)
}
