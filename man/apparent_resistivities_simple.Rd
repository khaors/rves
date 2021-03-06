% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ves_utilities.R
\name{apparent_resistivities_simple}
\alias{apparent_resistivities_simple}
\title{apparent_resistivities_simple}
\usage{
apparent_resistivities_simple(par, filt, spacing)
}
\arguments{
\item{par}{A numeric vector grouping real resistivities and layer thicknesses}

\item{filt}{A numeric vector with the values of the coefficients of the linear filter for an Schlumberger array}

\item{spacing}{A numeric vector with the values of the AB/2 spacing.}
}
\value{
A list with the following elements:This section includes all the functions
\itemize{
\item appres: A numeric vector with the values of the apparent resistivities
\item ab2: A numeric vector with the values of the electrode spacing
}
}
\description{
Function to calculate the apparent resistivities using a given Earth model.
}
\seealso{
Other utilities functions: \code{\link{apparent_resistivities}},
  \code{\link{convolution}}, \code{\link{krtrans}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
