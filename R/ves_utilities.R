#' @section utilities functions:
#'
#' This section includes all the functions required in other sections of this package.
#' The functions included in this section are:
#'
#' convolution, krtrans, apparent_resistivities, apparent_resistivities_simple
#'
#' @docType package
#' @name rves
NULL
#' @title
#' convolution
#' @description
#' Function to calculate the discrete convolution between two numeric vectors
#' @param x A numeric vector of length lx
#' @param y A numeric vector of length ly
#' @return
#' A numeric vector of length (lx+ly-1) with the results of the convolution between the x and y vectors
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family utilities functions
#' @export
convolution <- function(x, y) {
  lx <- length(x)
  ly <- length(y)
  cxy <- vector("numeric", length = lx+ly-1)
  for(i in 1:lx) {
    xi <- x[i]
    for(j in 1:ly) {
      ij <- i+j-1
      cxy[ij] <- cxy[ij] + xi * y[j]
    }
  }
  return(cxy)
}
#' @title
#' krtrans
#' @description
#' Function to calculate the resistivity transform of a given model of a layered earth
#' @param lambda The transform parameter
#' @param nlayers the number of layers in the model
#' @param thick A numeric vector with the thicknesses of the layers in the Earth model
#' @param rho A numeric vector with the real resistivities of the layers in the Earth model
#' @return
#' A numeric vector with the result of the resistivity transform
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family utilities functions
#' @export
krtrans <- function(lambda, nlayers, thick, rho){
  trans <- rho[nlayers]
  for(i in seq(nlayers-1,1,-1)) {
    tlt <- tanh(thick[i]*lambda)
    trans <- (trans + rho[i]*tlt)/(1+trans*tlt/rho[i])
  }
  return(trans)
}
#' @title
#' apparent_resistivities
#' @description
#' Function to calculate the apparent resistivities using a given Earth model.
#' @param rho A numeric vector with the values of the real resistivities of an Earth model
#' @param thick A numeric vector with the values of the thicknesses of the Earth model
#' @param filt A numeric vector with the values of the coefficients of the linear filter for an Schlumberger array
#' @param spacing A numeric vector with the values of the AB/2 spacing.
#' @return
#' A list with the following elements:This section includes all the functions
#' \itemize{
#' \item appres: A numeric vector with the values of the apparent resistivities
#' \item ab2: A numeric vector with the values of the electrode spacing
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family utilities functions
#' @importFrom pracma interp1
#' @export
apparent_resistivities <- function(rho, thick, filt, spacing) {
  ncoeff <- length(filt)
  ntpoints <- 40
  dy <- 0.48052648
  Tr <- vector("numeric", length = ntpoints)
  lambda <- vector("numeric", length = ntpoints)
  spa <- vector("numeric", length = ntpoints)
  nlayers <- length(rho)
  for (k in seq(1, ntpoints, 1)) {
    off <- -5.9+dy*k
    lambda[k] <- 1/exp(off)
    Tr[k] <- krtrans(lambda[k], nlayers, thick, rho)
    spa[k] <- exp(off)
  }
  #
  result <- convolution(Tr, filt)
  #print(result)
  #
  nappres <- length(result)-2*ncoeff
  appres <- vector("numeric", length = (nappres+1))
  ab2 <- vector("numeric", length = (nappres+1))
  for(i in 0:nappres){
    appres[i+1] <- result[i+ncoeff]
    off <- dy*(i+2) #reference  CALIB_A The 2 found by trial and error
    ab2[i+1] <- exp(off)
  }
  #
  if(!missing(spacing)){
    appres1 <- NULL
    #tryCatch(interp1(x = ab2, y = appres, xi = spacing),
    #         error = function(error){})
    appres1 <- interp1(x = ab2, y = appres, xi = spacing)
    ab2 <- spacing
    appres <- appres1
  }
  res <- list("ab2" = ab2, "appres" = appres)
  return(res)
}
#' @title
#' apparent_resistivities_simple
#' @description
#' Function to calculate the apparent resistivities using a given Earth model.
#' @param par A numeric vector grouping real resistivities and layer thicknesses
#' @param filt A numeric vector with the values of the coefficients of the linear filter for an Schlumberger array
#' @param spacing A numeric vector with the values of the AB/2 spacing.
#' @return
#' A list with the following elements:This section includes all the functions
#' \itemize{
#' \item appres: A numeric vector with the values of the apparent resistivities
#' \item ab2: A numeric vector with the values of the electrode spacing
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family utilities functions
#' @export
apparent_resistivities_simple <- function(par, filt, spacing){
  npar <- length(par)
  nres <- npar/2
  res <- apparent_resistivities(par[1:nres], par[(nres+1):npar], filt, spacing)
  return(res$appres)
}
