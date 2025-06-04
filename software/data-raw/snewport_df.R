library(surveillance)
library(lubridate)
library(ggplot2)

data("salmNewport")
salmNewportGermany <- aggregate(salmNewport, by = "unit")

snewport_df <- data.frame(
  date = epoch(salmNewportGermany),
  observed = salmNewportGermany@observed[, 1],
  id = "S. Newport Germany"
)

usethis::use_data(snewport_df, overwrite = TRUE)
