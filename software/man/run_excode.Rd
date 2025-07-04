% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_excode.R
\name{run_excode}
\alias{run_excode}
\alias{run_excode,sts,excodeModel-method}
\alias{run_excode,list,excodeModel-method}
\alias{run_excode,data.frame,excodeModel-method}
\title{Detect Excess Counts in Epidemiological Time Series}
\usage{
run_excode(
  surv_ts,
  excode_model,
  timepoints,
  learning_type = "unsupervised",
  maxIter = 100,
  verbose = FALSE,
  return_full_model = FALSE,
  past_timepoints_not_included = 26,
  time_units_back = 5
)

\S4method{run_excode}{sts,excodeModel}(
  surv_ts,
  excode_model,
  timepoints,
  learning_type = "unsupervised",
  maxIter = 100,
  verbose = FALSE,
  return_full_model = FALSE,
  past_timepoints_not_included = 26,
  time_units_back = 5
)

\S4method{run_excode}{list,excodeModel}(
  surv_ts,
  excode_model,
  timepoints,
  learning_type = "unsupervised",
  maxIter = 100,
  verbose = FALSE,
  return_full_model = FALSE,
  past_timepoints_not_included = 26,
  time_units_back = 5
)

\S4method{run_excode}{data.frame,excodeModel}(
  surv_ts,
  excode_model,
  timepoints,
  learning_type = "unsupervised",
  maxIter = 100,
  verbose = FALSE,
  return_full_model = FALSE,
  past_timepoints_not_included = 26,
  time_units_back = 5
)
}
\arguments{
\item{surv_ts}{A surveillance time series input. This can be one of the following:
\itemize{
  \item A \code{data.frame}, which can represent either a single or multiple time series.
        This data.frame contains the following columns:
        \itemize{
          \item \strong{date}: The date of aggregation.
          \item \strong{observed}: The observed number of cases.
          \item \strong{id}: The identifier of the time series (even if only one).
          \item \strong{offset} (optional): For example, the susceptible population.
          \item \strong{state} (optional): Indicates periods with special events such as outbreaks
            (e.g., 0 = normal, 1 = excess cases, NA = unknown).
        }
  \item A single \code{sts} object containing one timeseries.
  \item A named list of \code{sts} objects, each representing a separate time series.
}}

\item{excode_model}{An object of class \code{\linkS4class{excodeModel}} which specifies the model parameters and structure}

\item{timepoints}{An \code{integer} or a sequence of integers specifying
the time points for which excess count detection should be performed.}

\item{learning_type}{A \code{character} string indicating the type of learning.
Must be one of \code{c("unsupervised", "semisupervised", "supervised")}.}

\item{maxIter}{An \code{integer} specifying the maximum number of iterations
for the Expectation-Maximization (EM) algorithm. Defaults to \code{100}.}

\item{verbose}{A \code{logical} indicating whether progress should be printed
during execution. Defaults to \code{FALSE}.}

\item{return_full_model}{A \code{logical} indicating whether to return the full
model output. If \code{FALSE} (default), only results for the specified
\code{timepoints} are returned. If \code{TRUE}, the complete time series used
for model fitting is returned.}

\item{past_timepoints_not_included}{An \code{integer} specifying the number of past
time points to exclude from initialization. Defaults to \code{26}. If your data is daily then the past x days are excluded, for weekly data the past x weeks are excluded.}

\item{time_units_back}{An \code{integer} specifying how many past time units to include when fitting the model.
The default is \code{5}. For example, if \code{timepoints_per_unit = 52} (weekly data) was specified in excodeFormula,
then this corresponds to using the past 5 years of data. If \code{timepoints_per_unit = 7} (daily data and the units are weeks) was specified in excodeFormula then time_units_back equal to \code{5} corresponds to the past 5 weeks.}
}
\value{
Returns a fitted \code{\linkS4class{excodeModel}} object.
}
\description{
This function fits an `excodeModel` to epidemiological surveillance data
and returns the estimated parameters of an `excodeFamily`. It can handle
single or multiple time series.
}
\details{
This function currently accepts a single `sts` object, a `data.frame`, or a **named list of `sts` objects**.
However, a **an `sts` object matrix** such as the `measlesDE` dataset from the surveillance package is currently not supported and will result in an error.
}
\examples{

# Creating a Poisson harmonic model for the examples
excode_family_pois <- excodeFamily("Poisson")
excode_formula_har <- excodeFormula("Harmonic")
excode_har_pois <- excodeModel(
  excode_family_pois,
  excode_formula_har
)
# Example 1: Using data.frame as input time series
data(shadar_df)
result_shadar_har <- run_excode(shadar_df, excode_har_pois, 209:295)


# Example 2: Using an sts object (from 'surveillance' package) as input time series
library(surveillance)
data(stsNewport)
result_newport_har <- run_excode(stsNewport, excode_har_pois, 209:295)

# Example 3: Using a named list of two sts objects as input
stsShadar <- surveillance::sts(shadar_df$observed,
  epoch = shadar_df$date,
  state = shadar_df$state
)
named_list <- c("salmNewport" = stsNewport, "shadar" = stsShadar)
result_list <- run_excode(
  named_list,
  excode_har_pois, 290
)

}
\seealso{
\code{\linkS4class{excodeModel}}
}
