#' Fit Functional Relationship s = f(m)
#' Clause 8.5
#'
#' @param m Vector of mean levels
#' @param s Vector of standard deviations (s_r or s_R)
#' @param rel Relationship type: 'I' (s=bm), 'II' (s=a+bm), 'III' (s^2=a_v^2+(b_v m)^2), 'IV' (lg s = c + d lg m)
fit_precision_relationship <- function(m, s, rel = c("I", "II", "III", "IV"), max_iter = 10) {
  rel <- match.arg(rel)
  q <- length(m)

  if (rel == "I") {
    # Clause 8.5.2.4: Simplified formula b = sum(s/m) / q
    b <- sum(s / m) / q
    return(list(rel = "I", formula = "s = b * m", params = c(b = b), fitted = b * m))
  }

  if (rel == "II") {
    # Clause 8.5.2: Weighted iterative regression for s = a + b*m
    s_hat <- s
    for (iter in 1:max_iter) {
      W <- 1 / (s_hat^2)
      T1 <- sum(W)
      T2 <- sum(W * m)
      T3 <- sum(W * m^2)
      T4 <- sum(W * s)
      T5 <- sum(W * m * s)

      denom <- (T1 * T3 - T2^2)
      a <- (T3 * T4 - T2 * T5) / denom
      b <- (T1 * T5 - T2 * T4) / denom
      s_hat <- a + b * m
    }
    return(list(rel = "II", formula = "s = a + b * m", params = c(a = a, b = b), fitted = s_hat))
  }

  if (rel == "III") {
    # Clause 8.5.3: s^2 = a_v^2 + (b_v * m)^2
    s_hat <- s
    s_sq <- s^2
    m_sq <- m^2
    for (iter in 1:max_iter) {
      W_v <- 1 / (s_hat^4)
      Tv1 <- sum(W_v)
      Tv2 <- sum(W_v * m_sq)
      Tv3 <- sum(W_v * m_sq^2)
      Tv4 <- sum(W_v * s_sq)
      Tv5 <- sum(W_v * m_sq * s_sq)

      denom <- (Tv1 * Tv3 - Tv2^2)
      a_v_sq <- (Tv3 * Tv4 - Tv2 * Tv5) / denom
      b_v_sq <- (Tv1 * Tv5 - Tv2 * Tv4) / denom

      a_v <- sqrt(max(0, a_v_sq))
      b_v <- sqrt(max(0, b_v_sq))
      s_hat <- sqrt(a_v^2 + (b_v * m)^2)
    }
    return(list(rel = "III", formula = "s^2 = a_v^2 + (b_v * m)^2", params = c(a_v = a_v, b_v = b_v), fitted = s_hat))
  }

  if (rel == "IV") {
    # Clause 8.5.4: Unweighted regression lg s = c + d * lg m
    lg_m <- log10(m)
    lg_s <- log10(s)

    T1 <- sum(lg_m)
    T2 <- sum(lg_m^2)
    T3 <- sum(lg_s)
    T4 <- sum(lg_m * lg_s)

    denom <- q * T2 - T1^2
    c <- (T2 * T3 - T1 * T4) / denom
    d <- (q * T4 - T1 * T3) / denom

    C_const <- 10^c
    s_hat <- C_const * (m^d)
    return(list(rel = "IV", formula = "s = C * m^d", params = c(C = C_const, c = c, d = d), fitted = s_hat))
  }
}
