% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-excodeFormula-class.R
\docType{class}
\name{MultiState-class}
\alias{MultiState-class}
\title{This class is a container for the parameterization of a MultiState model.}
\description{
This class is a container for the parameterization of a MultiState model.
}
\section{Slots}{

\describe{
\item{\code{nStates}}{The number of states of the model.}

\item{\code{data}}{A data.frame containing variables which are used to model case counts. All variables in the data.frame will be used in the model. The data.frame has to have the same number of rows as time points in the time series.}

\item{\code{intercept}}{TRUE if the model should include an intercept, FALSE otherwise.}

\item{\code{timepoints_per_unit}}{Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).}

\item{\code{offset}}{TRUE if an offset should be included in the model.}

\item{\code{formula_bckg}}{A formula which models the 'normal' (background) states.}

\item{\code{formula}}{A formula which models which includes variable(s) to model 'excess' state(s).}
}}

