test_that("GLS gives similar results as the ISO method for functional relations", {

  # test data.frame
  precision_df <- data.frame(
    m   = c(3.94, 8.28, 14.18, 15.59, 20.41),
    s_R = c(0.171, 0.498, 0.400, 0.579, 0.637)
  )

  # tolerance for comparison
  tlrnc <- 1e-4

  # model 0
  mod0_iso <- with(precision_df, fit_precision_relationship(m, s_R, rel = "0"))
  mod0_gls <- fit_precision_gls(precision_df, rel = "0")
  expect_equal(mod0_iso$fitted, as.numeric(fitted(mod0_gls)),
               tolerance = tlrnc)

  # model 1
  mod1_iso <- with(precision_df, fit_precision_relationship(m, s_R, rel = "I"))
  mod1_gls <- fit_precision_gls(precision_df, rel = "I")
  expect_equal(mod1_iso$fitted, as.numeric(fitted(mod1_gls)),
               tolerance = tlrnc)

  # model 2
  mod2_iso <- with(precision_df, fit_precision_relationship(m, s_R, rel = "II"))
  mod2_gls <- fit_precision_gls(precision_df, rel = "II")
  expect_equal(mod2_iso$fitted, as.numeric(fitted(mod2_gls)),
               tolerance = tlrnc)

  # model 3
  mod3_iso <- with(precision_df, fit_precision_relationship(m, s_R, rel = "III"))
  mod3_gls <- fit_precision_gls(precision_df, rel = "III")
  expect_equal(mod3_iso$fitted, sqrt(as.numeric(fitted(mod3_gls))),
               tolerance = tlrnc)

  # model 4
  mod4_iso <- with(precision_df, fit_precision_relationship(m, s_R, rel = "IV"))
  mod4_gls <- fit_precision_gls(precision_df, rel = "IV")
  expect_equal(mod4_iso$fitted, 10^as.numeric(fitted(mod4_gls)),
               tolerance = tlrnc)

})
