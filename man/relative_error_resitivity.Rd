% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{relative_error_resitivity}
\alias{relative_error_resitivity}
\title{relative_error_resistivity}
\usage{
relative_error_resitivity(rho, thick, spacing, rho.measured)
}
\arguments{
\item{rho}{A numeric vector with the values of the real resistivity}

\item{thick}{A numeric vector with the values of the layer thicknesses}

\item{spacing}{A numeric vector with the values of the electrode spacing}

\item{rho.measured}{A numeric vector with the values of the apparent resistivity}
}
\value{
This function returns the value of the relative error.
}
\description{
Function to calculate the relative error between the calculated and measured resistivities.
}
\seealso{
Other calibration functions: \code{\link{calibrate_ilsqp}},
  \code{\link{calibrate_joint_nls}},
  \code{\link{calibrate_nls}},
  \code{\link{calibrate_seq_joint_nls}},
  \code{\link{calibrate_step_nls}},
  \code{\link{calibrate_step}},
  \code{\link{calibrate_svd}}, \code{\link{calibrate}},
  \code{\link{log_mnad_resistivity}},
  \code{\link{log_mxad_resistivity}},
  \code{\link{log_rss_resistivity}},
  \code{\link{mnad_resistivity}},
  \code{\link{mxad_resistivity}},
  \code{\link{rss_resistivity}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
