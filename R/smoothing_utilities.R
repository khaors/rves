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
#' @examples
#' library(ggplot2)
#' data("ves_data1")
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' #Smoothing with smooth.spline
#' res.ss <- smoothing_ves(sev1a, method = "smooth.spline")
#' sev1ss <- ves(id = "Sounding 1-Smoothing Spline", ab2 = res.ss$ab2,
#'               apprho = res.ss$apprho)
#' rho <- c(40,70,30, 20)
#' thick <- c(2,10,50,500)
#' par0 <- c(rho, thick)
#' res.ss.nls <- calibrate_nls(sev1ss, par0, iterations = 20, trace = FALSE)
#' sev1ss$rhopar <- res.ss.nls$rho
#' sev1ss$thickpar <- res.ss.nls$thickness
#' sev1ss$interpreted <- TRUE
#' p1 <- plot(sev1ss, type = "ves")
#' p1 <- p1 + geom_point(aes(x = ab2, y = apprho), data = ves_data1)
#' print(p1)
#' #Smoothing with kernel regression
#' res.kr <- smoothing_ves(sev1a, method = "kernel.regression", bw = 0.5)
#' sev1kr <- ves(id = "Sounding 1-Kernel Regression", ab2 = res.kr$ab2,
#'               apprho = res.kr$apprho)
#' res.kr.nls <- calibrate_nls(sev1kr, par0, iterations = 20, trace = FALSE)
#' sev1kr$rhopar <- res.kr.nls$rho
#' sev1kr$thickpar <- res.kr.nls$thickness
#' sev1kr$interpreted <- TRUE
#' p2 <- plot(sev1kr, type = "ves")
#' p2 <- p2 + geom_point(aes(x = ab2, y = apprho), data = ves_data1)
#' print(p2)
#' #Smoothing using wavelet thresholding
#' res.wv <- smoothing_ves(sev1a, method = "wavelet")
#' sev1wv <- ves(id = "Sounding 1-Wavelet Thresholding", ab2 = res.wv$ab2,
#'               apprho = res.wv$apprho)
#' res.wv.nls <- calibrate_nls(sev1wv, par0, iterations = 20, trace = FALSE)
#' sev1wv$rhopar <- res.wv.nls$rho
#' sev1wv$thickpar <- res.wv.nls$thickness
#' sev1wv$interpreted <- TRUE
#' p3 <- plot(sev1wv, type = "ves")
#' p3 <- p3 + geom_point(aes(x = ab2, y = apprho), data = ves_data1)
#' print(p3)
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
