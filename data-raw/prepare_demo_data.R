# prepare demo data for ringtest package

ringtest_demo <- read.csv(
  "data-raw/ringtest_data.csv"
)

usethis::use_data(
  ringtest_demo,
  overwrite = TRUE
)
