% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ves.R
\name{plot.ves}
\alias{plot.ves}
\title{plot.ves}
\usage{
\method{plot}{ves}(x, main = NULL, type = c("ves", "transformation",
  "diagnostic"), trans.type = c("direct", "scaling", "zohdy",
  "zohdy.smoothed"), ...)
}
\arguments{
\item{x}{A VES object}

\item{main}{Title of the plot}

\item{type}{A character string specifying the plot type. Currently only ves (measurements and
earth model), transformation (resistivity-depth transformation) and diagnostic are supported.}

\item{trans.type}{A character string specifying the transformation type. Only direct, scaling,
zohdy and zohdy.smoothed  are currently supported.}

\item{...}{Additional parameters to be passed to the plot function}
}
\description{
Function to plot the Vertical Electric Sounding data. This function can create two different
types of plots: resistivity and interpretation. The resistivity plot includes the
apparent resistivity vs electrode spacing. In the interpretation plot, the values of the
apparent resistivity vs electrode spacing and the thickness and real resisitivity of the
layers are included.
}
\examples{
library(grid)
library(gridExtra)
data(ves_data1)
ab2 <- ves_data1$ab2
apprho <- ves_data1$apprho
sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
# Conventional VES plot
plot(sev1a, type = "ves")
# Apply direct transformation to sev1a
plot(sev1a, type = "transformation", trans.type = "direct")
# Apply scaling transformation to sev1a
plot(sev1a, type = "transformation", trans.type = "scaling")
# Apply Zohdy transformation to sev1a
plot(sev1a, type = "transformation", trans.type = "zohdy")
# Apply Zohdy smoothed transformation to sev1a
plot(sev1a, type = "transformation", trans.type = "scaling")
# Diagnostic plot
rho <- c(40,70,30, 20)
thick <- c(2,10,50,500)
par0 <- c(rho, thick)
res.nls <- calibrate_nls(sev1a, par0, iterations = 30, ireport = 5, trace = FALSE)
sev1a$rhopar <- res.nls$rho
sev1a$thickpar <- res.nls$thickness
sev1a$interpreted <- TRUE
p1 <- plot(sev1a, type = "ves")
print(p1)
}
\seealso{
Other base functions: \code{\link{print.ves}},
  \code{\link{summary.ves}}, \code{\link{ves}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
