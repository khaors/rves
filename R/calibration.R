#' @section calibration functions:
#'
#' This section includes all the functions required for the estimation of the thicknesses and real resitivities from the VES data. The functions included here are:
#'
#' rss_resisitivity, log_rss_resistivity, mnad_resistivity, log_mnad_resistivity, mxad_resistivity, log_mxad_resistivity, calibrate
#'
#' @docType package
#' @name rves
NULL
#' @title
#' rss_resistivity
#' @description
#' Function to calculate the residual sum of squares between the measured and calculated apparent resistivity.
#' @param par A numeric vector with the values of the real resisitivity and the layer thickness
#' @param filter A numeric vector with the filter coefficients for a schlumberger array
#' @param apprho_measured A numeric vector with the values of the measured apparent resistivity
#' @param spacing A numeric vector with the values of AB/2.
#' @return
#' A numeric value with the residual sum of squares
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
rss_resistivity <-function(par, filter, apprho_measured, spacing){
  npar <- length(par)/2
  rho <- par[1:npar]
  thick <- par[(npar+1):length(par)]
  apprho_calculated <- apparent_resistivities(rho, thick, filter, spacing)
  #print(length(apprho_calculated$appres))
  res <- sqrt(mean((apprho_measured-apprho_calculated$appres)^2))
  return(res)
}
#' @title
#' log_rss_resistivity
#' @description
#' Function to calculate the log of the residual sum of squares between the measured and calculated apparent resistivity.
#' @param par A numeric vector with the values of the real resisitivity and the layer thickness
#' @param filter A numeric vector with the filter coefficients for a schlumberger array
#' @param apprho_measured A numeric vector with the values of the measured apparent resistivity
#' @param spacing A numeric vector with the values of AB/2.
#' @return
#' A numeric value with the log of the residual sum of squares
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
log_rss_resistivity <- function(par, filter, apprho_measured, spacing){
  npar <- length(par)/2
  rho <- par[1:npar]
  thick <- par[(npar+1):length(par)]
  apprho_calculated <- apparent_resistivities(rho, thick, filter, spacing)
  res <- sqrt(mean((log10(apprho_measured)-log10(apprho_calculated$appres))^2))
  return(res)
}
#' @title
#' mnad_resistivity
#' @description
#' Function to calculate the mean absolute deviation between the measured and calculated apparent resistivity.
#' @param par A numeric vector with the values of the real resisitivity and the layer thickness
#' @param filter A numeric vector with the filter coefficients for a schlumberger array
#' @param apprho_measured A numeric vector with the values of the measured apparent resistivity
#' @param spacing A numeric vector with the values of AB/2.
#' @return
#' A numeric value with the mean absolute deviation
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
mnad_resistivity <- function(par, filter, apprho_measured, spacing){
  npar <- length(par)/2
  rho <- par[1:npar]
  thick <- par[(npar+1):length(par)]
  apprho_calculated <- apparent_resistivities(rho, thick, filter, spacing)
  res <- mean(abs((apprho_measured-apprho_calculated$appres)))
  return(res)
}
#' @title
#' log_mnad_resistivity
#' @description
#' Function to calculate the log of the mean absolute deviation between the measured and calculated apparent resistivity.
#' @param par A numeric vector with the values of the real resisitivity and the layer thickness
#' @param filter A numeric vector with the filter coefficients for a schlumberger array
#' @param apprho_measured A numeric vector with the values of the measured apparent resistivity
#' @param spacing A numeric vector with the values of AB/2.
#' @return
#' A numeric value with the log of the mean absolute deviation
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
log_mnad_resistivity <- function(par, filter, apprho_measured, spacing){
  npar <- length(par)/2
  rho <- par[1:npar]
  thick <- par[(npar+1):length(par)]
  apprho_calculated <- apparent_resistivities(rho, thick, filter, spacing)
  res <- mean(abs((log10(apprho_measured)-log10(apprho_calculated$appres))))
}
#' @title
#' mxad_resistivity
#' @description
#' Function to calculate the maximum absolute deviation between the measured and calculated apparent resistivity.
#' @param par A numeric vector with the values of the real resisitivity and the layer thickness
#' @param filter A numeric vector with the filter coefficients for a schlumberger array
#' @param apprho_measured A numeric vector with the values of the measured apparent resistivity
#' @param spacing A numeric vector with the values of AB/2.
#' @return
#' A numeric value with the maximum absolute deviation
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
mxad_resistivity <- function(par, filter, apprho_measured, spacing){
  npar <- length(par)/2
  rho <- par[1:npar]
  thick <- par[(npar+1):length(par)]
  apprho_calculated <- apparent_resistivities(rho, thick, filter, spacing)
  res <- max(abs((apprho_measured-apprho_calculated$appres)))
  return(res)
}
#' @title
#' log_mxad_resistivity
#' @description
#' Function to calculate the log of the maximum absolute deviation between the measured and calculated apparent resistivity.
#' @param par A numeric vector with the values of the real resisitivity and the layer thickness
#' @param filter A numeric vector with the filter coefficients for a schlumberger array
#' @param apprho_measured A numeric vector with the values of the measured apparent resistivity
#' @param spacing A numeric vector with the values of AB/2.
#' @return
#' A numeric value with the log of the maximum absolute deviation
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
log_mxad_resistivity <- function(par, filter, apprho_measured, spacing){
  npar <- length(par)/2
  rho <- par[1:npar]
  thick <- par[(npar+1):length(par)]
  apprho_calculated <- apparent_resistivities(rho, thick, filter, spacing)
  res <- max(abs((log10(apprho_measured)-log10(apprho_calculated$appres))))
  return(res)
}
#' @title
#' calibrate
#' @description
#' Function to estimate the resistivities and thicnesses of the layers.
#' @param ves A VES object
#' @param opt.method A character string specifying the optimization method used to estimate the real resisitivity
#' @param obj.fn Objective function used in the parameter estimation
#' @param par0 A numeric vector with the Initial solution
#' @param lower A numeric vector of length equal to the number of layers with the min values of the parameter space
#' @param upper A numeric vector of length equal to the number of layers with the max values of the parameter space
#' @param control.par A list with the parameters of the optimization algorithm
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item a
#' \item b
#' \item c
#' }
#' @importFrom stats optim
#' @importFrom GenSA GenSA
#' @importFrom GA ga
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family calibration functions
calibrate <- function(ves, opt.method = c("L-BFGS-B", "SA", "GA"),
                      obj.fn = c('rss', 'mnad', 'mxad', 'log_rss', 'log_mnad', 'log_mxad'),
                      par0 = par0, lower = lower, upper = upper, control.par = NULL){
  if(class(ves) != "ves"){
    stop('ERROR: A ves object is required as input')
  }
  par.obj.fn <- NULL
  if(obj.fn == 'rss'){
    par.obj.fn <- function(par, filter, apprho_measured, spacing){
      rss_resistivity(par, filter, apprho_measured, spacing)
    }
  }
  else if(obj.fn == 'mnad'){
    par.obj.fn <- function(par, filter, apprho_measured, spacing){
      mnad_resistivity(par, filter, apprho_measured, spacing)
    }
  }
  else if(obj.fn == 'mxad'){
    par.obj.fn <- function(par, filter, apprho_measured, spacing){
      mxad_resistivity(par, filter, apprho_measured, spacing)
    }
  }
  else if(obj.fn == 'log_rss'){
    par.obj.fn <- function(par, filter, apprho_measured, spacing){
      log_rss_resistivity(par, filter, apprho_measured, spacing)
    }
  }
  else if(obj.fn == "log_mnad"){
    par.obj.fn <- function(par, filter, apprho_measured, spacing){
      log_mnad_resistivity(par, filter, apprho_measured, spacing)
    }
  }
  else if(obj.fn == 'log_mxad'){
    par.obj.fn <- function(par, filter, apprho_measured, spacing){
      log_mxad_resistivity(par, filter, apprho_measured, spacing)
    }
  }
  res1 <- NULL
  if(opt.method == "L-BFGS-B"){
    res <- optim(par0,  fn = par.obj.fn, filter = as.matrix(rves::filt$V1),
                 apprho_measured = ves$appres, spacing = ves$ab2,
                 method = 'L-BFGS-B', lower = lower, upper = upper)
    res1 <- list(par = res$par, value = res$value)
  }
  else if(opt.method == "SA"){
    if(is.null(control.par)){
      control.par <- list(maxit = 100, threshold.stop = 0.1, nb.stop.improvement = 50,
                          smooth = FALSE, simple.function = TRUE, verbose = TRUE)
    }
    res <- GenSA(par0, fn = par.obj.fn, lower = lower, upper = upper,
                     control = control.par, filter = as.matrix(rves::filt$V1),
                     apprho_measured = ves$appres, spacing = ves$ab2)
    res1 <- list(par = res$par, value = res$value)
  }
  else if(opt.method == "GA"){
    if(is.null(control.par)){
      control.par <- list(method = "L-BFGS-B", poptim = 0.05, pressel = 0.5,
                          control = list(fnscale = -1, maxit = 100))
    }
    res <- ga(type = "real-valued",
                 fitness =  function(par, filter, apprho_measured, spacing) {
                   -par.obj.fn(par, filter, apprho_measured, spacing)
                 },
                 filter = as.matrix(rves::filt$V1), apprho_measured = ves$appres,
                 spacing = ves$ab2, min = lower, max = upper,
                 popSize = 300, maxiter = 500, pcrossover = 0.85, pmutation = .2,
                 run = 50, maxFitness = -.05,
                 optim = TRUE, optimArgs = control.par)
    print(res)
    res1 <- list(par = res@solution, value = res@bestvalue)
  }
  return(res1)
}
#' @title
#' calibrate_nls
#' @description
#' Function to estimate the real resistivities and thicknesses using the nonlinear least-squares approach
#' proposed by Roy(1999).
#' @param ves A ves object
#' @param par0 A numeric vector with the values of the layer resistivities and thicknesses
#' @param iterations Number of iterations
#' @param ireport Number of iterations to report results on the console
#' @return
#' A list with the following entries:
#' \itemize{
#' \item par: A numeric vector with the values of the layer resistivities and thicknesses
#' \item value: The value or the RSS (Residual sum of squares)
#' \item rel.erro: The value of the relative error (in percentage)
#' \item cal.error: A matrix with the RSS and relative error at each iteration
#' \item hessian: First order approximation of the Hessian matrix. This is calculate using the Jacobian matrix.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#' Roy, I. An efficient non-linear least-squares 1D inversion scheme for resistivity
#' and IP sounding data. 1999. Geophysical Prospecting, 47, 4, 527-550.
#' @family calibration functions
#' @importFrom numDeriv jacobian
#' @importFrom pracma mod
#' @export
calibrate_nls <- function(ves, par0, iterations = 100, ireport = 10){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  #
  cal.error <- matrix(0.0, nrow = iterations, ncol = 2)
  spacing <- ves$ab2
  measured <- ves$appres
  npar <- length(par0)
  nparh <- npar/2
  niter <- iterations
  iter <- 0
  cpar <- par0
  cpar1 <- as.matrix(cpar, npar, 1)
  mu_max <- 20
  mu <- mu_max
  beta <- 2
  while(iter < niter){
    # Data difference
    cres <- apparent_resistivities_simple(cpar, rves::filt$V1, spacing)
    d <- as.matrix((cres-measured), length(measured), 1)
    # Calculate jacobian matrix
    J <- jacobian(apparent_resistivities_simple, cpar, filt=rves::filt$V1, spacing=spacing)
    # calculate parameter correction
    x <- solve(t(J)%*%J+mu*pracma::eye(npar),t(J)%*%d)#/(iter+1)
    # calculate error vector
    err <- J%*%x-d
    # Calculate y
    y <- solve(t(J)%*%J+(mu)*pracma::eye(npar),t(J)%*%err)
    E1 <- norm(y)
    #
    cpar1 <- cpar1 - x
    cpar <- as.numeric(cpar1)
    cres <- apparent_resistivities_simple(cpar, rves::filt$V1, spacing)
    current.error <- mean((cres-measured)**2)
    current.error1 <- 100*mean(abs(measured-cres)/measured)
    cal.error[iter+1,1] <- mean((cres-measured)**2)
    cal.error[iter+1,2] <- 100*mean(abs(measured-cres)/measured)
    if(iter != -1){
      mu2 <- mu/2
      x2 <- solve(t(J)%*%J+mu2*pracma::eye(npar),t(J)%*%d)
      err2 <- J%*%x2-d
      y2 <- solve(t(J)%*%J+mu2*pracma::eye(npar),t(J)%*%err2)
      E2 <- norm(y2)
      p <- sqrt(E1/E2)
      if(p > 1){
        q <- exp(-beta*mu2)-p*exp(-beta*mu)
        r <- abs((1-p)/q)
        new_mu <- log(r)/beta
      }
      else {
        new_mu <- mu2
      }
      if((mu-1)>1e-4){
        mu <- new_mu
      }
    }
    if(mod(iter,ireport) == 0 | iter == (niter-1)){
      cat("iteration, RSS, Rel Error = ", c(iter,current.error,current.error1), "\n")
    }
    iter <- iter+1
  }
  J <- jacobian(apparent_resistivities_simple, cpar, filt=rves::filt$V1, spacing=spacing)
  res <- list(par = cpar, value = current.error, rel.error = current.error1,
              cal.error = cal.error, rho = cpar[1:nparh],
              thickness = cpar[(nparh+1):npar], hessian = t(J)%*%J)
  return(res)
}
