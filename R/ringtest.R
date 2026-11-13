#' ringtest
#'
#' @description
#' Calculate general mean and variance components for each level (ISO clause 8.6.11)
#'
#' @param input_data Two types of `input_data` are accecpted. The data must be cleaned first and outliers removed. Option (1) A data.frame in long format of individual replicate measurements with columns `level`, `lab`, `value`; the REML-procedure will be used to calculate the ring-test results. Option (2) A data.frame with aggregated cell-statistics (using the function `calc_cell_stats`) with columns  `level`, `lab`, `n`, `mean`, `sd`; the classic ISO procedure will be used to calculate ring-test results.
#'
#' @return A data.frame with columns `level`, `p` (number of labs), `m_hat` (general mean), `s_r` (repeatability SD), `s_L` (between-lab SD), `s_R` (reproducibility SD)
#'
#' @export
ringtest <- function(input_data) {

  # check type of input_data: individual or aggregated, error otherwise
  cols   <- colnames(input_data)
  c_aggr <- c("level", "lab", "n", "mean", "sd")
  c_ind  <- c("level", "lab", "value")
  type   <- ifelse(
    sum(c_ind %in% cols) == 3, "individual",
    ifelse(sum(c_aggr %in% cols) == 5, "aggregated",
           "format not recognized")
  )
  stopifnot("format of 'input_data' not recognized" =
              type %in% c("individual", "aggregated"))

  levels <- unique(input_data$level)
  precision_tbl <- do.call(rbind, lapply(levels, function(lvl) {
    input_data_sub <- input_data[input_data$level == lvl, ]
    res <- switch(
      type,
      individual = calc_variances_reml(input_data_sub),
      aggregated = calc_variances_iso(input_data_sub)
    )
    res$level <- lvl
    res
  }))
}
