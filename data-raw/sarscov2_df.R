library(readr)
library(lubridate)

# Source: https://robert-koch-institut.github.io/SARS-CoV-2-Infektionen_in_Deutschland/

# TODO: get this csv file from Bene as it is publicly available and store it inside the data-raw folder
tab <- read_csv("C:/Users/zacherb/Downloads/Aktuell_Deutschland_SarsCov2_Infektionen.csv")
# then replace with this
tab <- read_csv("data-raw/Aktuell_Deutschland_SarsCov2_Infektionen.csv")

tab_count <- tab %>%
  filter(IdLandkreis %in% 11008) %>%
  group_by(Meldedatum) %>%
  count()

tab_neukölln <- data.frame(Meldedatum = seq(min(tab_count$Meldedatum), as_date("2020-07-31"), by = 1))
sarscov2_df <- left_join(tab_neukölln, tab_count) %>%
  mutate(
    id = "SARS-CoV-2 Berlin-Neukölln",
    n = ifelse(is.na(n), 0, n)
  ) %>%
  rename(
    date = Meldedatum,
    observed = n
  )

usethis::use_data(sarscov2_df, overwrite = TRUE)
