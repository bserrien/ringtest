#' fit_precision_gls
#'
#' @description
#' Fit Functional Relationship s_R = f(m) (ISO Clause 8.5)
#'
#' @importFrom nlme gls varPower
#'
#' @param data data.frame containing variables `m` and `s_R` (outcome of function precision)
#' @param rel Relationship type:
#' '0' (s = constant) (8.6.13.2)
#' 'I' (s = b * m),
#' 'II' (s = a + b * m),
#' 'III' (s^2 = a_v^2 + (b_v * m)^2),
#' 'IV' (lg s = c + d * lg m)
#'
#' @export
fit_precision_gls <- function(data, rel) {

  stopifnot(
    "argument 'rel' should be one from '0', 'I', 'II', 'III', 'IV'" =
      rel %in% c("0", "I", "II", "III", "IV")
  )

  # fit chosen model with generalized least squares
  if (rel == "0") {
    model <- gls(
      s_R ~ 1,
      method = "REML",
      data = data
    )
  } else if (rel == "I") {
    model <- gls(
      s_R ~ -1 + m,
      weights = varPower(fixed = 1, form = ~ m),
      method = "REML",
      data = data
    )
  } else if (rel == "II") {
    model <- gls(
      s_R ~ 1 + m,
      weights = varPower(fixed = 1, form = ~ fitted(.)),
      method = "REML",
      data = data
    )
  } else if (rel == "III") {
    data$m_sq   <- data$m^2
    data$s_R_sq <- data$s_R^2
    model <- gls(
      s_R_sq ~ 1 + m_sq,
      weights = varPower(fixed = 1, form = ~ fitted(.)),
      method = "REML",
      data = data
    )
  } else if (rel == "IV") {
    model <- gls(
      log10(s_R) ~ log10(m),
      method = "REML",
      data = data
    )
  }
  return(model)
}
