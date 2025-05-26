library(surveillance)
library(lubridate)
library(ggplot2)

data(shadar)

shadar_df <- data.frame(
  date = as_date("2001-01-01") + (shadar$week) * 7,
  observed = shadar$observed[, 1],
  id = "S. Hadar Germany"
)

usethis::use_data(shadar_df, overwrite = TRUE)
