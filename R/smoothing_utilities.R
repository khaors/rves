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
#' @param bw Bandwidth
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
#' @importFrom stats smooth.spline approx
#' @importFrom locpol locCuadSmootherC EpaK
#' @importFrom wavethresh wd wr threshold
#' @importFrom pracma interp1
#' @export
smoothing_ves <- function(ves, method = c("smooth.spline", "kernel.regression",
                                          "wavelet"),
                          bw = 0.1, ...){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  ab2 <- ves$ab2
  appres <- ves$appres
  ndat <- length(ab2)
  #print(appres)
  #
  res <- NULL
  if(method == "smooth.spline"){
    tmp <- smooth.spline(log10(ab2), log10(appres))
    ab2.s <- 10^tmp$x
    rho.s <- 10^tmp$y
    if(sum(is.na(ab2.s)) >= 1 | sum(is.na(rho.s)) >= 1){
      warning("The smoothing method has numerical instabilities. No smoothing applied.")
      ab2.s <- ab2
      rho.s <- appres
    }
    res <- list(ab2 = ab2.s, apprho = rho.s)
  }
  else if(method == "kernel.regression"){
    lab2 <- log10(ab2)
    lapprho <- log10(appres)
    tmp <- locCuadSmootherC(lab2, lapprho, lab2, bw = bw, EpaK)
    ab2.s <- 10^tmp$x
    rho.s <- 10^tmp$beta0
    if(sum(is.na(tmp$x)) >= 1 | sum(is.na(tmp$beta0)) >= 1){
      warning("The smoothing method has numerical instabilities. No smoothing applied.")
      ab2.s <- ab2
      rho.s <- appres
    }
    res <- list(ab2 = ab2.s, apprho = rho.s)
  }
  else if(method == "wavelet"){
    current.power <- floor(log(ndat)/log(2))+1
    ab2.out <- logseq(min(log10(ab2)), 0.99*max(log10(ab2)), 2**current.power)
    apprho.approx <- approx(log10(ab2),log10(appres),
                            xout = log10(ab2.out), method = 'linear')
    pos.valid <- is.na(apprho.approx$y)
    #if(sum(pos.valid) > 0)
    waveletwmap <- wd(apprho.approx$y, family="DaubLeAsymm", filter.number=10)
    softthreshwmap <- threshold(waveletwmap, type="soft", policy="universal")
    hardthreshwmap <- threshold(waveletwmap, type="hard", policy="universal")
    s.soft <- wr(softthreshwmap)
    s.hard <- wr(hardthreshwmap)
    xout <- c(1.01*ab2[1],ab2[2:(ndat-1)],0.99*max(ab2.out))
    soft <- interp1(ab2.out, s.soft, xi = xout)
    hard <- interp1(ab2.out, s.hard, xi = xout)
    if(sum(is.na(soft)) >= 1){
      warning("The smoothing method has numerical instabilities. No smoothing applied.")
      soft <- log10(appres)
    }
    res <- list(ab2 = ab2, apprho = 10^soft)
  }
  return(res)
}
#' @title
#' logseq
#' @description
#' Function to generate sequences with logarithmic increments
#' @param from Initial point of the sequence
#' @param to Final point of the sequence
#' @param n Number of points
#' @return
#' This function returns a numeric vector with the generated sequence.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
logseq <- function( from = 1, to = 1, n) {
  exp(log(10)*seq(from, to, length.out = n))
}
