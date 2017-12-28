#' @section calibration functions:
#'
#' This section includes all the functions required for the estimation of the thicknesses and real resitivities from the VES data. The functions included here are:
#'
#' rss_resisitivity, log_rss_resistivity, mnad_resistivity, log_mnad_resistivity, mxad_resistivity, log_mxad_resistivity, calibrate, calibrate_nls, calibrate_svd
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
#' relative_error_resistivity
#' @description
#' Function to calculate the relative error between the calculated and measured resistivities.
#' @param rho A numeric vector with the values of the real resistivity
#' @param thick A numeric vector with the values of the layer thicknesses
#' @param spacing A numeric vector with the values of the electrode spacing
#' @param rho.measured A numeric vector with the values of the apparent resistivity
#' @return
#' This function returns the value of the relative error.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
relative_error_resitivity <- function(rho, thick, spacing, rho.measured){
  rho.calc <- apparent_resistivities(rho, thick, rves::filt$V1, spacing)
  rel.err <- mean(100*abs(rho.calc$appres - rho.measured)/rho.measured)
  return(rel.err)
}
#' @title
#' calibrate
#' @description
#' Function to estimate the resistivities and thicnesses of the layers.
#' @param ves A VES object
#' @param opt.method A character string specifying the optimization method used to estimate the real resisitivity.
#' Currently the following methods are supported:
#' \itemize{
#' \item L-BFGS-B: Limited-memory modification of the BFGS quasi-Newton method (optim package)
#' \item SA: Simulated Annealing (GenSA package)
#' \item GA: Genetic Algorithms (GA package)
#' \item PSO: Particle Swarm Optimization (pso package)
#' \item DE: Differential Evoluation (DEoptim package)
#' }
#' @param obj.fn Objective function used in the parameter estimation
#' @param par0 A numeric vector with the Initial solution
#' @param lower A numeric vector of length equal to the number of layers with the min values of the parameter space
#' @param upper A numeric vector of length equal to the number of layers with the max values of the parameter space
#' @param control.par A list with the parameters of the optimization algorithm
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item par: vector with all parameters
#' \item value: value of the objective function
#' \item rho: Numeric vector with the real resistivities
#' \item thickness: Numeric vector with the layer thicknesses
#' \item rel.error: The value of the relative error
#' }
#' @importFrom stats optim
#' @importFrom GenSA GenSA
#' @importFrom GA ga
#' @importFrom pso psoptim
#' @importFrom DEoptim DEoptim DEoptim.control
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @family calibration functions
calibrate <- function(ves, opt.method = c("L-BFGS-B", "SA", "GA", "PSO", "DE"),
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
  nparh <- length(par0)/2
  npar <- length((par0))
  if(opt.method == "L-BFGS-B"){
    res <- optim(par0,  fn = par.obj.fn, filter = as.matrix(rves::filt$V1),
                 apprho_measured = ves$appres, spacing = ves$ab2,
                 method = 'L-BFGS-B', lower = lower, upper = upper)
    rho = res$par[1:nparh]
    thick = res$par[(nparh+1):npar]
    spacing <- ves$ab2
    rel.err <- relative_error_resitivity(rho, thick, spacing, ves$appres)
    res1 <- list(par = res$par, value = res$value, rho = rho,
                 thickness = thick, rel.err = rel.err)
  }
  else if(opt.method == "SA"){
    if(is.null(control.par)){
      control.par <- list(maxit = 100, threshold.stop = 0.1, nb.stop.improvement = 50,
                          smooth = FALSE, simple.function = TRUE, verbose = TRUE)
    }
    res <- GenSA(par0, fn = par.obj.fn, lower = lower, upper = upper,
                     control = control.par, filter = as.matrix(rves::filt$V1),
                     apprho_measured = ves$appres, spacing = ves$ab2)
    rho = res$par[1:nparh]
    thick = res$par[(nparh+1):npar]
    spacing <- ves$ab2
    rel.err <- relative_error_resitivity(rho, thick, spacing, ves$appres)
    res1 <- list(par = res$par, value = res$value, rho = rho,
                 thickness = thick, rel.err = rel.err)
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
                 popSize = 300, maxiter = 100, pcrossover = 0.85, pmutation = .2,
                 run = 50, maxFitness = -.05,
                 optim = TRUE, optimArgs = control.par)
    rho <- res@solution[1:nparh]
    thick <- res@solution[(nparh+1):npar]
    spacing <- ves$ab2
    rel.err <- relative_error_resitivity(rho, thick, spacing, ves$appres)
    res1 <- list(par = res@solution, value = abs(res@fitnessValue),
                 rho = res@solution[1:nparh],
                 thickness = res@solution[(nparh+1):npar], rel.err = rel.err)
  }
  else if(opt.method == "PSO"){
    if(is.null(control.par)){
      control.par <- list(trace = 0, maxit.stagnate = 30, maxit = 300)
    }
    res.pso <- psoptim(par = par0, fn = par.obj.fn,
                       filter = as.matrix(rves::filt$V1),
                       apprho_measured = ves$appres,
                       spacing = ves$ab2,
                       lower = lower,
                       upper= upper, control = control.par)
    rho <- res.pso$par[1:nparh]
    thick <- res.pso$par[(nparh+1):npar]
    spacing <- ves$ab2
    rel.err <- relative_error_resitivity(rho, thick, spacing, ves$appres)
    res1 <- list(par = res.pso$par, value = res.pso$value, rho = rho,
                 thickness = thick, rel.err = rel.err)
  }
  else if(opt.method == "DE"){
    if(is.null(control.par)){
      control.par <- DEoptim.control(strategy=1, NP=100, itermax=300, CR=0.9,
                                     F=0.8, trace=0, storepopfrom=1, storepopfreq=1)
    }
    res.de <- suppressWarnings(DEoptim(par.obj.fn,
                                       lower = lower,
                                       upper = upper,
                                       control = control.par,
                                       filter = as.matrix(rves::filt$V1),
                                       apprho_measured = ves$appres,
                                       spacing = ves$ab2))
    rho <- unname(res.de$optim$bestmem[1:nparh])
    thick <- unname(res.de$optim$bestmem[(nparh+1):npar])
    spacing <- ves$ab2
    rel.err <- relative_error_resitivity(rho, thick, spacing, ves$appres)
    res1 <- list(par = res.de$optim$bestmem, value = res.de$optim$bestval,
                 rho = rho, thickness = thick, rel.err = rel.err)
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
#' \item rel.error: The value of the relative error (in percentage)
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
#' @title
#' calibrate_svd
#' @description
#' Function to estimate the layer parameters using the approach proposed by Meju (1992)
#' @param ves A vertical electrical sounding object
#' @param par0 A numeric vector with the initial values of resistivities and thicknesses
#' @param iterations An integer specifying the maximum number of iterations
#' @param ireport An integer specifying the report interval
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item par: A numeric vector with the values of the layer resistivities and thicknesses
#' \item value: The value or the RSS (Residual sum of squares)
#' \item rel.error: The value of the relative error (in percentage)
#' \item cal.error: A matrix with the RSS and relative error at each iteration
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family Calibration functions
#' @references Meju, M. An effective ridge regression procedure for resistivity data inversion.
#' Computers and Geosciences, 18(2-3), 99-119, 1992.
#' @importFrom numDeriv jacobian
#' @importFrom pracma inv
#' @export
calibrate_svd <- function(ves, par0, iterations = 100, ireport = 10){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  cal.error <- matrix(0.0, nrow = iterations, ncol = 2)
  spacing <- ves$ab2
  measured <- ves$appres
  npar <- length(par0)
  nparh <- npar/2
  if(nparh <= 2){
    stop('ERROR: the number of layers used is small. Increase the number of layers')
  }
  niter <- iterations
  iter <- 0
  cpar <- par0
  cpar1 <- as.matrix(cpar, npar, 1)
  ks <- seq(from = 1,to = 10,by = 1)
  qk <- vector('numeric', length = 10)
  current.error1 <- 100
  while(iter < niter){
    # Data difference
    cres <- apparent_resistivities_simple(cpar, rves::filt$V1, spacing)
    rel.err <- mean(100*abs(cres-measured)/measured)
    if(rel.err < 10.0 & iter < 1){
      break
    }
    d <- as.matrix((cres-measured), length(measured), 1)
    # Calculate jacobian matrix
    J <- jacobian(apparent_resistivities_simple, cpar, filt=rves::filt$V1, spacing=spacing)
    # Calculate SVD of Jacobian Matrix
    Jsvd <- svd(J)
    # Get the eigenvalues of the Jacobian matrix
    qvalues <- Jsvd$d
    #
    if(min(qvalues)<1.e-5){
      p <- which.min(qvalues)
      qvalues[p] <- 1e-5
    }
    # Define the search range for the optimun q
    ql <- 10*max(qvalues)
    qs <- 0.1*min(qvalues)
    qk <- ((100*qs-ql)+(ql-qs)*ks^2)/99.0
    beta <- qk^2
    # calculate parameter correction
    #x <- solve(t(J)%*%J+beta[ik]*pracma::eye(npar),t(J)%*%d)
    rel.err1 <- vector('numeric', length = 11)
    rel.err1[1] <- rel.err
    for(ii in 1:10){
      for(ik in (11-ii):1){
        J <- jacobian(apparent_resistivities_simple, cpar, filt=rves::filt$V1,
                      spacing=spacing)
        Jsvd <- svd(J)
        corrected.eigenvalues <- diag(Jsvd$d + beta[ik])
        x <- Jsvd$v%*%pracma::inv(corrected.eigenvalues)%*%t(Jsvd$u)%*%d
        oldpar <- cpar
        cpar <- cpar - x
        cres1 <- apparent_resistivities_simple(cpar, rves::filt$V1, spacing)
        rel.err1[ii+1] <- mean(100*abs(cres1-measured)/measured)
        #print(rel.err1[ii+1])
        if(rel.err1[ii] > rel.err || rel.err1[ii+1] > rel.err1[ii]){
          cpar <- oldpar
          break
        }
        else{
          cpar <- cpar
          current.error1 <- rel.err1[ii+1]
          current.error <- sum((cres1-measured)^2)
        }
      }
    }
    if(mod(iter, ireport) == 0 | iter == (niter-1)){
      cat("iteration, RSS, Rel Error = ", c(iter,current.error,current.error1), "\n")
    }
    iter <- iter + 1
  }
  current.error <- 100.0
  res <- list(par = cpar, value = current.error, rel.error = current.error1,
              cal.error = cal.error, rho = cpar[1:nparh],
              thickness = cpar[(nparh+1):npar])
  return(res)
}
#' @title
#' calibration_ilsqp
#' @description
#' Function to estimate the layer parameters (real resistivities and thicknesses) from
#' a VES measurements using the ILPSQ (Iterative Least-Squares Procedure) with Singular
#' Value Decomposition proposed by Muiuane and Pedersen (1999).
#' @param ves A VES object
#' @param iterations An integer specifying the maximum number of iterations
#' @param ireport An integer specifying the report interval
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item par: A numeric vector with the values of the layer resistivities and thicknesses
#' \item value: The value or the RSS (Residual sum of squares)
#' \item rel.error: The value of the relative error (in percentage)
#' \item cal.error: A matrix with the RSS and relative error at each iteration
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family Calibration functions
#' @references
#' Muiuane, E. & Pedersen, L. Automatic 1D interpretation of DC resistivity sounding
#' data. Journal of Applied Geophysics, 42, 1, 35-45, 1999.
calibration_ilsqp <- function(ves, iterations = 100, ireport = 10){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  # Define initial model: two layers
  par0 <- c( mean(ves$appres), mean(ves$appres), ves$ab2[1]/4, 500)
  #print(par0)
  res.2layer <- calibrate_nls(ves, par0, iterations = iterations, ireport = ireport)
  res.current.layer <- res.2layer
  res.old.layer <- res.2layer
  current.error <- 100
  max.error <- res.2layer$rel.error
  while(current.error > max.error){
    nparh <- length(res.current.layer$par)/2
    new.rho <- vector("numeric", length = (nparh + 1))
    new.thick <- vector("numeric", length = (nparh + 1))
    #
    new.rho[1:(nparh-1)] <- res.current.layer$rho[1:(nparh-1)]
    new.rho[(nparh+1)] <- res.current.layer$rho[nparh]
    new.rho[nparh] <- res.current.layer$rho[(nparh-1)]
    #
    new.thick[1:(nparh-1)] <- res.current.layer$thickness[1:(nparh-1)]
    new.thick[(nparh+1)] <- res.current.layer$thickness[nparh]
    new.thick[nparh] <- res.current.layer$thickness[(nparh-1)]
    #
    new.par0 <- c(new.rho, new.thick)
    #print(new.par0)
    res.old.layer <- res.current.layer
    res.current.layer <- calibrate_nls(sev1, par0 = new.par0,
                                       iterations = iterations,
                                       ireport = ireport)
    #print(res.current.layer)
    max.error <- current.error
    current.error <- res.current.layer$rel.error
    cat("Current Error= ", current.error, "\n")
    #stop('ERROR')
  }
  res.current.layer <- calibrate_nls(ves, par0 = res.current.layer$par,
                                     iterations = iterations,
                                     ireport = ireport)
  return(res.current.layer)
}
