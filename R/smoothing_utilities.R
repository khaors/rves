#' @section smoothing functions:
#'
#' This section includes all the functions required for the smoothing of the VES data. The functions included here are:
#'
#' smoothing_ves
#'
#' @docType package
#' @name rves
NULL
#' @title
#' smoothing_ves
#' @description
#' Function to smooth the measurements of the apparent resistivity of a VES
#' @param ves A Vertical Electric Sounding objet
#' @param method A character string specifying the smoothing method. Currently the following methods are implemented:
#' \itemize{
#' \item smooth.spline
#' \item kernel.regression
#' }
#' @param ... Additional parameters required for the smoothing method
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item ab2: Values of the electrode separation
#' \item apprho: Values of the smoothed apparent resistivity
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family smoothing functions
#' @importFrom stats smooth.spline
#' @export
smoothing_ves <- function(ves, method = c("smooth.spline", "kernel.regression"), ...){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  ab2 <- ves$ab2
  appres <- ves$appres
  #
  res <- NULL
  if(method == "smooth.spline"){
    res <- smooth.spline(log10(ab2), log10(appres))
  }
  else if(method == "kernel.regression"){

  }
  return(res)
}
