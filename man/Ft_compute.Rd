% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ft_compute.R
\name{Ft_compute}
\alias{Ft_compute}
\title{function to compute the cumulative risk function Ft from incidence rates}
\usage{
Ft_compute(rates, rates_denom = 1e+05, smo = TRUE, smo_df = 15)
}
\arguments{
\item{rates}{Populational incidence rates by 1-year from 0y to 120y}

\item{rates_denom}{Number "at risk" for the incidence rates (eg 100,000)}

\item{smo}{Logical, smooth the Ft function or not ? Without smoothing the Ft function can reflect the age classes.}

\item{smo_df}{degree of freedom for the smooth function}
}
\value{
the cumulative risk function Ft from incidence rate data
}
\description{
function to compute the cumulative risk function Ft from incidence rates
}
\examples{
#to be completed
}
