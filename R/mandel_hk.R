#' calc_mandel_hk
#'
#' @description
#' Calculate Mandel's h and k statistics (ISO 8.3.2)
#'
#' @param cell_stats outcome of function `calc_cell_stats`
#'
#' @export
calc_mandel_hk <- function(cell_stats) {
  cell_stats$h <- NA_real_
  cell_stats$k <- NA_real_

  for (lvl in unique(cell_stats$level)) {
    idx <- which(cell_stats$level == lvl & !is.na(cell_stats$mean))
    p_j <- length(idx)

    # Formula 8.3.2.2: Mandel's h
    y_bar_ij          <- cell_stats$mean[idx]
    y_bar_j           <- mean(y_bar_ij) # grand mean at level j
    s_means           <- stats::sd(y_bar_ij)
    cell_stats$h[idx] <- (y_bar_ij - y_bar_j) / s_means

    # Formula 8.3.2.3: Mandel's k
    s_ij              <- cell_stats$sd[idx]
    s_ij_valid        <- s_ij[!is.na(s_ij)]
    pooled_s          <- sqrt(sum(s_ij_valid^2) / p_j)
    cell_stats$k[idx] <- cell_stats$sd[idx] / pooled_s
  }

  return(cell_stats)
}
