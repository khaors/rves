#' @section Transformation functions:
#'
#' This section includes functions used to calculate the resistivity-depth transformations, The results of these functions can be used for exploratory purposes or to define an initial solution of the automatic inversion procedures. .
#'
#' the functions included in this section are:
#'
#' transform_direct, transform_scaling, transform_zohdy, transform_smoothed_zohdy
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
#' @title
#' transform_zohdy
#' @description
#' Function used to obtain interpreted depths and resistivities from shifted electrode
#' spacings and adjusted apparent resistivities, respectively. This method does not
#' require optimization.
#' @param ves A VES object
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item depth: A numeric vector with the layer depths
#' \item real.res: A numeric vector with the layer resistivities
#' \item thickness: A numeric vector with the layer thicknesses
#' \item rho: A numeric vector with the layer resistivities (defined to be compatible
#' with calibration results)
#' \item rms.depth: A numeric value with the rms obtained for the depth conversion
#' \item rms.rho: A numeric value with the rms obtained for the resistivities.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family transformation functions
#' @references
#' Zohdy, A. A new method for the automatic interpretation of Schlumberger and Wenner
#' sounding curves. Geopysics, 1989, 54(2), 245-253.
#' @export
transform_zohdy <- function(ves){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  # a. ) Assume that the digitized electrode spacings are equal to
  # the depths and that the apparent resistivities are equal to the true
  # resistivities at those depths (Figure
  #
  obs.rho <-  ves$appres
  true.rho <- ves$appres
  current.apprho <- obs.rho
  thick <- ves$ab2
  spacing <- ves$ab2
  #
  # (b) Compute a theoretical sounding curve for this multilayer model by convolution
  #
  old.rho <- current.apprho
  current.apprho <- apparent_resistivities_simple(par = c(true.rho, thick),
                                                  rves::filt$V1,
                                                  spacing)
  #
  # (c) Compute the root-mean-square (rms) percent from the equation
  #
  current.rms <- sqrt(sum(((obs.rho-current.apprho)/obs.rho)^2)/length(obs.rho))
  #
  # (d) Multiply all the depths by 0.9 to decrease all the layer
  # depths by 10 percent (small arbitrary amount).
  #
  old.rms <- 1.5*current.rms
  while(old.rms > current.rms){
    old.thick <- thick
    thick <- 0.9*thick
    current.apprho <- apparent_resistivities_simple(par = c(true.rho, thick),
                                                    rves::filt$V1,
                                                    spacing)
    old.rms <- current.rms
    current.rms <- sqrt(sum(((obs.rho-current.apprho)/obs.rho)^2)/length(obs.rho))
  }
  rms.depth <- current.rms
  #
  # Resistivity correction sptep
  #
  #(a) At each digitized electrode spacing on the observed and calculated curves,
  #if the computed apparent resistivity, at the jth spacing, is less (or greater) than
  #the corresponding ob- served apparent resistivity, the corresponding true resistivity
  #of the jth layer should be increased (or decreased) so that the calculated apparent
  #resistivity will rise (or fall) and approach the observed resistivity (Figure 4d).
  #
  old.rms <- 1.5*current.rms
  thick <- old.thick
  while(old.rms > current.rms){
    current.apprho <- apparent_resistivities_simple(par = c(true.rho, thick),
                                                    rves::filt$V1,
                                                    spacing)
    true.rho <- true.rho*(obs.rho/current.apprho)
    old.rms <- current.rms
    current.rms <- sqrt(sum(((obs.rho-current.apprho)/obs.rho)^2)/length(obs.rho))
  }
  rms.rho <- current.rms
  depth <- cumsum(thick)
  real.res <- true.rho
  res <- list(depth = depth, real.res = real.res, thickness = thick, rho = real.res,
              rms.depth = rms.depth, rms.rho = rms.rho)
  return(res)
}
#' @title
#' transform_smoothed_zohdy
#' @description
#' Function used to improve the original Zohdy's method.
#' @param ves A VES object
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item depth: A numeric vector with the layer depths
#' \item real.res: A numeric vector with the layer resistivities
#' \item thickness: A numeric vector with the layer thicknesses
#' \item rho: A numeric vector with the layer resistivities (defined to be compatible
#' with calibration results)
#' \item rms.depth: A numeric value with the rms obtained for the depth conversion
#' \item rms.rho: A numeric value with the rms obtained for the resistivities.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family transformation functions
#' @references
#' Loke, M, \& Barker, R. Improvements to the Zohdy method for the inversion of
#' resistivity sounding and pseudosection data. Computers \& Geosciences, 21(2),
#' 321-332, 1995.
#' @export
transform_smoothed_zohdy <- function(ves){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  # a. ) Assume that the digitized electrode spacings are equal to
  # the depths and that the apparent resistivities are equal to the true
  # resistivities at those depths (Figure
  #
  obs.rho <-  ves$appres
  true.rho <- ves$appres
  current.apprho <- obs.rho
  thick <- ves$ab2
  spacing <- ves$ab2
  #
  # (b) Compute a theoretical sounding curve for this multilayer model by convolution
  #
  old.rho <- current.apprho
  current.apprho <- apparent_resistivities_simple(par = c(true.rho, thick),
                                                  rves::filt$V1,
                                                  spacing)
  #
  # (c) Compute the root-mean-square (rms) percent from the equation
  #
  current.rms <- sqrt(sum(((obs.rho-current.apprho)/obs.rho)^2)/length(obs.rho))
  #
  # (d) Multiply all the depths by 0.9 to decrease all the layer
  # depths by 10 percent (small arbitrary amount).
  #
  old.rms <- 1.5*current.rms
  while(old.rms > current.rms){
    old.rms <- current.rms
    old.thick <- thick
    thick <- 0.9*thick
    current.apprho <- apparent_resistivities_simple(par = c(true.rho, thick),
                                                    rves::filt$V1,
                                                    spacing)
    current.rms <- sqrt(sum(((obs.rho-current.apprho)/obs.rho)^2)/length(obs.rho))
  }
  rms.depth <- current.rms
  #
  # Resistivity correction sptep
  #
  #(a) At each digitized electrode spacing on the observed and calculated curves,
  #if the computed apparent resistivity, at the jth spacing, is less (or greater) than
  #the corresponding ob- served apparent resistivity, the corresponding true resistivity
  #of the jth layer should be increased (or decreased) so that the calculated apparent
  #resistivity will rise (or fall) and approach the observed resistivity (Figure 4d).
  #
  old.rms <- 1.5*current.rms
  thick <- old.thick
  Ccoeffs <- c(.25, .5, .25)
  while(old.rms > current.rms){
    current.apprho <- apparent_resistivities_simple(par = c(true.rho, thick),
                                                    rves::filt$V1,
                                                    spacing)
    current.e <- obs.rho/current.apprho
    smoothed.e <- current.e
    for(i in 2:(length(smoothed.e)-1)){
      smoothed.e[i] <- sum(Ccoeffs * c(current.e[i-1], current.e[i], current.e[i+1]))
    }
    true.rho <- true.rho*smoothed.e
    old.rms <- current.rms
    current.rms <- sqrt(sum(((obs.rho-current.apprho)/obs.rho)^2)/length(obs.rho))
  }
  #
  rms.rho <- current.rms
  depth <- cumsum(thick)
  real.res <- true.rho
  res <- list(depth = depth, real.res = real.res, thickness = thick, rho = real.res,
              rms.depth = rms.depth, rms.rho = rms.rho)
  #
  return(res)
}
