% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1-excodeFamily-class.R
\name{excodeFamily}
\alias{excodeFamily}
\title{Create a family of probability distributions for excess count detection.}
\usage{
excodeFamily(name, nb_size = NA, shared_nb_size = FALSE)
}
\arguments{
\item{name}{Name of the probability distribution that should be used. Either "Poisson" or "NegBinom".}

\item{nb_size}{Size parameter of the Negative Binomial distribution. Only relevant if 'MultiState' model is used with a Negative Binomial dsitribution.}

\item{shared_nb_size}{Logical indicating whether nb_size parameter is shared across multiple time series. Only relevant if name is "NegBinom".}
}
\value{
An \code{\linkS4class{excodeFamily}} object.
}
\description{
Create a family of probability distributions for excess count detection.
}
\examples{

excode_family_pois <- excodeFamily("Poisson")
excode_family_pois

}
\seealso{
\code{\linkS4class{excodeFamily}}
}
