#' @section Transformation functions:
#'
#' This section includes functions used to calculate the resistivity-depth transformations, The results of these functions can be used for exploratory purposes or to define an initial solution of the automatic inversion procedures. .
#'
#' the functions included in this section are:
#'
#' transform_direct, transform_scaling
#'
#' @docType package
#' @name rves
NULL
#' @title
#' transform_direct
#' @description
#' Function to transform the apparent resistivity into effective resistivities which are assumed closer to the real resistivities using an exponential equation.
#' @param ves A ves object.
#' @param k Empirical constant.
#' @param alpha Decay exponent (pseudo-porosity factor).
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item depth: Depth
#' \item real.res: Real resisitivity
#' }
#' @details
#' The conversion from apparent resitivities to effective resistivities is achieved using the following expression:
#' \deqn{
#' \rho_{eff}=k\rho_{app}\exp{(-(1-\alpha))}
#' }
#' where:
#' \itemize{
#' \item \eqn{\rho_{eff}} is the effective resistivity which is assumed to be closer to the true resistivity.
#' \item \eqn{\rho_{app}} is the apparent resistivity measured in the field.
#' \item \eqn{k} is a empirical constant defined to improve the fit between the apparent and effective resistivities. This constant has a value of 2.3.
#' \item \eqn{\alpha} is the exponent in the exponential relationship. Its value is of 0.15.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family transformation functions
#' @references
#' Meju, M.  Simple effective resistivity-depth transformations
#' for infield or real-time data processing. Computers and Geosciences,
#' 21(8), 985-992, 1995.
#' @export
transform_direct <- function(ves, k = 2.3, alpha = 0.15){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  real.res <- k*ves$appres*exp(-(1-alpha))
  depth <- ves$ab2/2.3
  res <- list(depth = depth, real.res = real.res)
  return(res)
}
#' @title
#' transform_scaling
#' @description
#' Function to
#' @param ves A VES object
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item depth: Depth
#' \item real.res: Real resisitivity
#' }
#' @details
#' The conversion from apparent to effective resistivities using the scaling approach is based on the use of the following expression:
#' \deqn{
#' \rho_{eff}^{i}=S \cdot \rho_{app}^{i}
#' }
#' where $S$ is the so called scaling factor defined as:
#' \deqn{
#' S = \frac{\rho_{app}^{i}}{ \rho_{app}^{i-1} }
#' }
#' and \eqn{\rho_{app}^{i}} and \eqn{\rho_{app}^{i-1}} are the apparent resistivities of the \eqn{i-}th and \eqn{(i-1)}-th layers respectively.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family transformation functions
#' @references
#' Meju, M.  Simple effective resistivity-depth transformations
#' for infield or real-time data processing. Computers and Geosciences,
#' 21(8), 985-992, 1995.
#' @export
transform_scaling <- function(ves){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  depth <- ves$ab2/2.3
  real.res <- vector("numeric", length = length(ves$ab2))
  real.res[1] <- ves$appres[1]
  for(ires in 2:length(ves$ab2)){
    S <- ves$appres[ires]/ves$appres[ires-1]
    real.res[ires] <- S*ves$appres[ires]
  }
  res <- list(depth = depth, real.res = real.res)
  return(res)
}
#
