#' Calculate Cell Means and Within-Cell Standard Deviations
#'
#' @param data A data.frame with columns: 'lab', 'level', 'value'
calc_cell_stats <- function(data) {
  labs <- unique(data$lab)
  levels <- unique(data$level)

  res <- expand.grid(lab = labs, level = levels, stringsAsFactors = FALSE)
  res$n <- 0
  res$mean <- NA_real_
  res$sd <- NA_real_

  for (i in seq_len(nrow(res))) {
    sub <- data[data$lab == res$lab[i] & data$level == res$level[i], "value"]
    sub <- sub[!is.na(sub)]
    res$n[i] <- length(sub)
    if (length(sub) > 0) {
      res$mean[i] <- mean(sub) # Formula (3)/(4) per cell
    }
    if (length(sub) > 1) {
      res$sd[i] <- stats::sd(sub) # Formula (1)
    }
  }
  return(res)
}

#' Calculate Mandel's h and k Statistics
#' Clause 8.3.2
calc_mandel_hk <- function(cell_stats) {
  cell_stats$h <- NA_real_
  cell_stats$k <- NA_real_

  for (lvl in unique(cell_stats$level)) {
    idx <- which(cell_stats$level == lvl & !is.na(cell_stats$mean))
    p_j <- length(idx)

    y_bar_ij <- cell_stats$mean[idx]
    y_bar_j <- mean(y_bar_ij) # Grand mean
    s_means <- stats::sd(y_bar_ij)

    # Formula 8.3.2.2: Mandel's h
    cell_stats$h[idx] <- (y_bar_ij - y_bar_j) / s_means

    # Formula 8.3.2.3: Mandel's k
    s_ij <- cell_stats$sd[idx]
    s_ij_valid <- s_ij[!is.na(s_ij)]
    pooled_s <- sqrt(sum(s_ij_valid^2) / p_j)
    cell_stats$k[idx] <- cell_stats$sd[idx] / pooled_s
  }
  return(cell_stats)
}
