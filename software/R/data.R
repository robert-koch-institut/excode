#' German all-cause mortality data
#'
#' A dataset containing the number of weekly deaths in Germany reported to the German Federal Statistical Office \url{https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Sterbefaelle-Lebenserwartung/sterbefallzahlen.html}. The variables are as follows:
#'
#' @format A tibble with 418 rows and 7 variables:
#' \describe{
#'   \item{date}{Week of death}
#'   \item{observed}{Number of deaths}
#'   \item{timepoint}{A numeric column counting sequential weeks from 0 to total number of weeks, representing time progression for time-trend modeling}
#'   \item{sin52}{Computed sinus component used for harmonic modeling}
#'   \item{cos52}{Computed cosine component used for harmonic modeling}
#'   \item{bckg_week}{A logical (`TRUE`/`FALSE`) variable marking background mortality weeks, used for estimating model parameters.}
#' }
"mort_df_germany"


#' SARS-CoV-2 infections in Berlin-Neukölln (Germany)
#'
#' A dataset containing the daily number of reported SARS-CoV-2 cases from March-July 2020 in Berlin-Neukölln \url{https://robert-koch-institut.github.io/SARS-CoV-2-Infektionen_in_Deutschland/}. The dataset was downloaded on 2024-10-31. The variables are as follows:
#'
#'
#' @format A data.frame with 154 rows and 3 variables:
#' \describe{
#'   \item{date}{Day of data collection}
#'   \item{observed}{Number of cases}
#' }
"sarscov2_df"

#' Salmonella Hadar cases in Germany 2001-2006
#'
#' A dataset containing the weekly number of reported Salmonella Hadar cases from January 2001 to August 2006 in Germany. The variables are as follows:
#'
#'
#' @format A data.frame with 295 rows and 3 variables:
#' \describe{
#'   \item{date}{Week of death}
#'   \item{observed}{Number of cases}
#' }
"shadar_df"
