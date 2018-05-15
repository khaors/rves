#' @section calibration functions:
#'
#' This section includes all the functions required for the estimation of the thicknesses and real resitivities from the VES data. The functions included here are:
#'
#' rss_resisitivity, log_rss_resistivity, mnad_resistivity, log_mnad_resistivity, mxad_resistivity, log_mxad_resistivity, calibrate, calibrate_nls, calibrate_svd,
#' calibrate_step_nls, calibrate_step
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
#' @param trace A logical flag to indicate if the optimization information must be traced
#' @return
#' A list with the following entries:
#' \itemize{
#' \item par: A numeric vector with the values of the layer resistivities and thicknesses
#' \item value: The value or the RSS (Residual sum of squares)
#' \item rel.error: The value of the relative error (in percentage)
#' \item cal.error: A matrix with the RSS and relative error at each iteration
#' \item residual: A vector with the residuals (original scale)
#' \item var.residual: The variance of the residuals (log transformed scale)
#' \item hessian: First order approximation of the Hessian matrix. This is calculated using the Jacobian matrix.
#' \item cov.matrix: Covariance matrix of the estimated parameters calculated as the inverse of the
#' hessian matrix.
#' \item corr.matrix: Correlation matrix of the estimated parameters.
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @references
#' Roy, I. An efficient non-linear least-squares 1D inversion scheme for resistivity
#' and IP sounding data. 1999. Geophysical Prospecting, 47, 4, 527-550.
#' @family calibration functions
#' @importFrom numDeriv jacobian
#' @importFrom pracma mod
#' @importFrom stats var
#' @export
#' @examples
#' # Example 1
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' # Let's fit a four layer model
#' rho <- c(40, 70, 30, 20)
#' thick <- c(2, 10, 50, 500)
#' par0 <- c(rho, thick)
#' res.nls1 <- calibrate_nls(sev1a, par0, iterations = 30, ireport = 5)
#' sev1a$rhopar <- res.nls1$rho
#' sev1a$thickpar <- res.nls1$thickness
#' sev1a$interpreted <- TRUE
#' plot(sev1a, type = "ves")
#' # Example 2
#' data(ves_data2)
#' ab2 <- ves_data2$ab2
#' apprho <- ves_data2$apprho
#' sev2a <- ves(id = "TD76", ab2 = ab2, apprho = apprho)
#' rho <- c(20, 50, 100, 10)
#' thick <- c(10, 20, 500, 100)
#' par0 <- c(rho, thick)
#' res.nls2 <- calibrate_nls(sev2a, par0, iterations = 30, ireport = 5)
#' sev2a$rhopar <- res.nls2$rho
#' sev2a$thickpar <- res.nls2$thickness
#' sev2a$interpreted <- TRUE
#' plot(sev2a, type = "ves")
calibrate_nls <- function(ves, par0, iterations = 30, ireport = 10, trace = TRUE){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  #
  f_model <- function(x, spacing){
    par <- 10^x
    res <- apparent_resistivities_simple(par, rves::filt$V1, spacing)
    #pos.valid <- !is.na(res) & res > 0.0
    #res <- res[pos.valid]
    #tryCatch(log10(res), finally = print(res))
    return(log10(res))
  }
  #
  cal.error <- matrix(0.0, nrow = iterations, ncol = 2)
  spacing <- ves$ab2
  measured <- log10(ves$appres)
  npar <- length(par0)
  nparh <- npar/2
  niter <- iterations
  iter <- 0
  cpar <- log10(par0)
  cpar1 <- as.matrix(cpar, npar, 1)
  mu_max <- 20
  mu <- mu_max
  beta <- 2
  old.par <- NULL
  old.error <- 100
  current.error1 <- 100
  while(iter < niter){
    # Data difference
    cres <- f_model(cpar, spacing)
    pos.invalid <- is.na(cres)
    if(sum(pos.invalid) >= 1){
      current.error1 <- old.error
      cpar <- old.par
      break
    }
    d <- as.matrix((cres-measured), length(measured), 1)
    # Calculate jacobian matrix
    J <- jacobian(f_model, x = cpar, spacing = spacing)
    # calculate parameter correction
    x <- solve(t(J)%*%J+mu*pracma::eye(npar),t(J)%*%d)#/(iter+1)
    # calculate error vector
    err <- J%*%x-d
    # Calculate y
    y <- solve(t(J)%*%J+(mu)*pracma::eye(npar),t(J)%*%err)
    E1 <- norm(y)
    #
    old.par <- cpar
    cpar1 <- cpar1 - x
    cpar <- as.numeric(cpar1)
    cres <- f_model(cpar, spacing)
    pos.invalid <- is.na(cres)
    if(sum(pos.invalid) >= 1){
      current.error1 <- old.error
      cpar <- old.par
      break
    }
    current.error <- mean((cres-measured)**2)
    old.error <- current.error1
    current.error1 <- 100*sqrt(mean(((measured-cres)/measured)^2))
    cal.error[iter+1,1] <- current.error
    cal.error[iter+1,2] <- current.error1
    if(current.error1 > old.error){
      if(trace){
        cat("Final step ", iter, ", Error= ", old.error, "\n")
      }
      current.error1 <- old.error
      cpar <- old.par
      break
    }
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
    if(trace){
      if(mod(iter,ireport) == 0 | iter == (niter-1)){
        cat("iteration, RSS, Rel Error = ", c(iter,current.error,current.error1), "\n")
      }
    }
    iter <- iter+1
  }
  J <- jacobian(f_model, x = cpar, spacing = spacing)
  cres <- f_model(cpar, spacing)
  residual.res <- (10^measured -10^cres)
  var.residual <- var(measured-cres) #var(residual.res)
  cov.matrix <- pracma::inv(t(J)%*%J + 1e-10*pracma::eye(length(cpar)))*var.residual
  nr <- nrow(cov.matrix)
  corr.matrix <- matrix(0.0, nrow = nr, ncol = nr)
  sd.par <- sqrt(diag(cov.matrix))
  for(i in 1:nr){
    for(j in i:nr){
      corr.matrix[i,j] <- cov.matrix[i,j]/(sd.par[i]*sd.par[j])
      corr.matrix[j,i] <- corr.matrix[i,j]
    }
  }
  res <- list(par = 10^cpar, value = current.error,
              rel.error = current.error1,
              cal.error = cal.error, rho = 10^cpar[1:nparh],
              thickness = 10^cpar[(nparh+1):npar],
              hessian = t(J)%*%J,
              cov.matrix = cov.matrix,
              residual = residual.res,
              var.residual = var.residual,
              corr.matrix = corr.matrix)
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
#' @param trace A logical flag to indicate if the optimization information must be traced
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
#' @family calibration functions
#' @references Meju, M. An effective ridge regression procedure for resistivity data inversion.
#' Computers and Geosciences, 18(2-3), 99-119, 1992.
#' @importFrom numDeriv jacobian
#' @importFrom pracma inv
#' @export
calibrate_svd <- function(ves, par0, iterations = 100, ireport = 10, trace = TRUE){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  #
  f_model <- function(x, spacing){
    par <- 10^x
    res <- apparent_resistivities_simple(par, rves::filt$V1, spacing)
    return(log10(res))
  }
  #
  cal.error <- matrix(0.0, nrow = iterations, ncol = 2)
  spacing <- ves$ab2
  measured <- log10(ves$appres)
  npar <- length(par0)
  nparh <- npar/2
  if(nparh <= 2){
    stop('ERROR: the number of layers used is small. Increase the number of layers')
  }
  niter <- iterations
  iter <- 0
  cpar <- log10(par0)
  cpar1 <- as.matrix(cpar, npar, 1)
  ks <- seq(from = 1,to = 10,by = 1)
  qk <- vector('numeric', length = 10)
  current.error1 <- 100
  while(iter < niter){
    # Data difference
    #cres <- apparent_resistivities_simple(cpar, rves::filt$V1, spacing)
    cres <- f_model(cpar, spacing)
    rel.err <- mean(100*abs(10^cres-10^measured)/10^measured)
    if(rel.err < 1.0 & iter < 1){
      break
    }
    d <- as.matrix((cres-measured), length(measured), 1)
    # Calculate jacobian matrix
    #J <- jacobian(apparent_resistivities_simple, cpar, filt=rves::filt$V1, spacing=spacing)
    J <- jacobian(f_model, x = cpar, spacing = spacing)
    # Calculate SVD of Jacobian Matrix
    Jsvd <- svd(J)
    # Get the eigenvalues of the Jacobian matrix
    qvalues <- Jsvd$d
    # Define the search range for the optimun q
    ql <- 10*max(qvalues)
    qs <- 0.1*min(qvalues)
    qk <- ((100*qs-ql)+(ql-qs)*ks^2)/99.0
    beta <- qk^2
    #print(beta)
    # calculate parameter correction
    #x <- solve(t(J)%*%J+beta[ik]*pracma::eye(npar),t(J)%*%d)
    rel.err1 <- vector('numeric', length = 11)
    rel.err1[1] <- rel.err
    current.error1 <- 100
    current.error <- 100
    for(ii in 1:10){
      for(ik in (11-ii):1){
        corrected.eigenvalues <- diag(qvalues + beta[ik])
        if(determinant(corrected.eigenvalues, logarithm = TRUE)$modulus < 1e-5){
          break
        }
        x <- Jsvd$v%*%pracma::inv(corrected.eigenvalues)%*%t(Jsvd$u)%*%d
        oldpar <- cpar
        cpar <- cpar - x
        cres1 <- f_model(cpar, spacing)
        rel.err1[ii+1] <- mean(100*abs(10^cres1-10^measured)/10^measured)
        if(rel.err1[ii] > rel.err || rel.err1[ii+1] > rel.err1[ii]){
          cpar <- oldpar
          current.error1 <- rel.err1[ii+1]
          current.error <- sum((cres1-measured)^2)
          break
        }
        else{
          cpar <- cpar
          current.error1 <- rel.err1[ii+1]
          current.error <- sum((10^cres1-10^measured)^2)
        }
      }
    }
    if(trace){
      if(mod(iter, ireport) == 0 | iter == (niter-1)){
        cat("iteration, RSS, Rel Error = ", c(iter,current.error,current.error1), "\n")
      }
    }
    iter <- iter + 1
  }
  current.error <- 100.0
  res <- list(par = 10^cpar, value = current.error, rel.error = current.error1,
              cal.error = cal.error, rho = 10^cpar[1:nparh],
              thickness = 10^cpar[(nparh+1):npar])
  return(res)
}
#' @title
#' calibrate_ilsqp
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
#' @family calibration functions
#' @export
#' @references
#' Muiuane, E. & Pedersen, L. Automatic 1D interpretation of DC resistivity sounding
#' data. Journal of Applied Geophysics, 42, 1, 35-45, 1999.
calibrate_ilsqp <- function(ves, iterations = 30, ireport = 10){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  rho.mn <- min(ves$appres)/3
  rho.mx <- 3*max(ves$appres)
  # Define initial model: two layers
  par0 <- c( mean(ves$appres), mean(ves$appres), ves$ab2[1]/4, 500)
  #print(par0)
  res.2layer <- calibrate_nls(ves, par0, iterations = iterations, ireport = ireport)
  #res.2layer <- calibrate(ves, opt.method = "L-BFGS-B", obj.fn = "rss",
  #                         par0 = par0, lower = c(rep(rho.mn,2), rep(5,2)),
  #                         upper = c(rep(rho.mx,2),rep(501,2)))
  res.current.layer <- res.2layer
  res.old.layer <- res.2layer
  current.error <- 100
  max.error <- res.2layer$rel.err
  nparh <- 0
  best.model <- NULL
  min.error <- 100
  #print(max.errorsev1$rhopar<-tmp1$rho
  for(ilayer in 1:19){
    old.error <- max.error
    #while(current.error > max.error){
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
      #
      res.current.layer <- calibrate_nls(ves, par0 = new.par0, iterations = iterations,
                                         ireport = ireport)
      #res.current.layer <- calibrate(ves, opt.method = "L-BFGS-B", obj.fn = "log_rss",
      #                               par0 = new.par0,
      #                               lower = c(rep(rho.mn, nparh), rep(5,nparh)),
      #                               upper = c(rep(rho.mx,nparh), rep(501,nparh)))
      #print(res.current.layer)
      if(res.current.layer$rel.err < min.error){
        min.error <- res.current.layer$rel.err
        best.model <- res.current.layer
        cat("Current best model, Nr Layers= ", as.character(ilayer + 1), "\n")
      }
      current.error <- res.current.layer$rel.err
      cat("Current Error= ", current.error, "\n")
      #stop('ERROR')
    #}
  }
  #
  #res.current.layer <- calibrate(ves, opt.method = "L-BFGS-B", obj.fn = "rss",
  #                               par0 = res.current.layer$par,
  #                               lower = c(rep(rho.mn, nparh), rep(5,nparh)),
  #                               upper = c(rep(rho.mx, nparh), rep(501,nparh)))
  res.current.layer <- calibrate_nls(ves, par0 = best.model$par,
                                     iterations = iterations,
                                     ireport = ireport)
  return(res.current.layer)
}
#' @title
#' calibrate_step_nls
#' @description
#' Stepwise estimation of true resistivities and thicknesses of a VES.
#' @param ves A VES object
#' @param iterations An integer specifying the maximum number of iterations
#' @param ireport An integer specifying the report interval
#' @param max.layers An integer specifying the maximum number of layers to include in the
#' sequential estimation.
#' @param select.measure A character string specifying the criteria used to select the best
#' model from the set of models tried during the sequential estimation. Available options
#' are:
#' \itemize{
#' \item rss: (default) sum of residual squares
#' \item aic: Akaike Information Criteria
#' \item bic: Bayesian Information Criteria
#' }
#' @param trace A logical flag to indicate if the optimization information must be traced
#' @return
#' This function returns a list with the following entries:
#' \itemize{
#' \item par: A numeric vector with the values of the layer resistivities and thicknesses
#' \item value: The value or the RSS (ResiduSpatial Data Analysis in Ecology and Agricultureal sum of squares)
#' \item rel.error: The value of the relative error (in percentage)
#' \item cal.error: A matrix with the RSS and relative error at each iteration
#' }
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family calibration functions
#' @export
calibrate_step_nls <- function(ves, iterations = 30, ireport = 10,
                              max.layers = 10, select.measure = "rss",
                              trace = TRUE){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  #
  current.res <- NULL
  current.par <- NULL
  best.res <- NULL
  current.err <- 0.1
  current.aic <- 1e10
  max.aic <- 1e12
  max.bic <- 1e12
  max.err <- 100
  all.measures <- matrix(0.0, nrow = (max.layers-1), ncol = 3)
  depth <- ves$ab2/2.3
  thick <- diff(depth)
  n <- length(depth)
  for(ilay in 2:max.layers){
    cat(paste0("Current model: ", as.character(ilay), " Layers \n"))
    pos <- vector("numeric", length = ilay+1)
    if(ilay <= length(thick)){
     pos <- round(seq(1, length(thick), length.out = (ilay+1)))
     #pos <- round(seq(1, length(thick), length.out = (ilay)))
    }
    else{
     pos <- c(seq(1, length(thick)), rep(length(thick), (ilay-length(thick)+1) ))
    }
    current.par <- c(rep(mean(ves$appres), ilay), thick[pos])
    current.res <- calibrate_nls(ves, par0 = current.par, iterations = iterations,
                                 ireport = ireport, trace = trace)
    current.err <- current.res$rel.error
    current.lik <- -(n/2)*log(2*pi)-(n/2)*log(current.res$value)
    current.aic <- 2*(2*ilay-1)-2*current.lik
    current.bic <- log(n)*(2*ilay-1)-2*current.lik
    all.measures[(ilay-1),] <- c(current.err, current.aic, current.bic)
    if(select.measure == "rss"){
      if(current.err < max.err){
        best.res <- current.res
        max.err <- current.err
      }
    }
    else if(select.measure == "aic"){
      if(current.aic < max.aic){
        best.res <- current.res
        max.aic <- current.aic
      }
    }
    else if(select.measure == "bic"){
      if(current.bic < max.bic){
        best.res <- current.res
        max.bic <- current.bic
      }
    }
  }
  best.res$all.measures <- all.measures
  res <- best.res
  return(best.res)
}
#' @title
#' calibrate_step
#' @description
#' Stepwise VES inversion using different optimization algorithms
#' @param ves A VES object
#' @param opt.method A character string specifying the optimization method to be used
#' @param max.layers An integer with the maximum number of layers to be testsed
#' @param lower A numeric value
#' @param upper A numeric value
#' @param select.measure  A character string specifying the criteria used to select the best
#' model from the set of models tried during the sequential estimation. Available options
#' are:
#' \itemize{
#' \item rss: (default) sum of residual squares
#' \item aic: Akaike Information Criteria
#' \item bic; Bayesian Information Criteria
#' }
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
#' @family calibration functions
#' @export
calibrate_step <- function(ves, opt.method, max.layers = 10, lower = 1, upper = 500,
                           select.measure = "rss"){
  if(class(ves) != "ves"){
    stop('ERROR: A VES object is required as input')
  }
  current.res <- NULL
  current.par <- NULL
  best.res <- NULL
  current.err <- 0.1
  current.aic <- 1e10
  max.aic <- 1e12
  max.bic <- 1e12
  max.err <- 100
  all.measures <- matrix(0.0, nrow = (max.layers-1), ncol = 3)
  depth <- ves$ab2/2.3
  thick <- diff(depth)
  n <- length(depth)
  cat(paste0("Optimization Method: ", opt.method, "\n"))
  for(ilay in 2:max.layers){
    if(ilay <= length(thick)){
      pos <- round(seq(1, length(thick), length.out = ilay))
    }
    else{
      pos <- c(seq(1, length(thick)), rep(length(thick), (ilay-length(thick)+1) ))
    }
    current.par <- c(rep(mean(ves$appres), ilay), thick[pos])
    cat(paste0("Current model: ", as.character(ilay), " Layers \n"))
    #print(ilay)
    #print(pos)
    #print(current.par)
    lowerlim <- rep(rep(lower, ilay), 2)
    upperlim <- rep(rep(upper, ilay), 2)
    current.res <- calibrate(ves, par0 = current.par,
                             opt.method = opt.method,
                             obj.fn = "log_rss",
                             lower = lowerlim,
                             upper = upperlim)
    current.err <- current.res$rel.err
    current.lik <- -(n/2)*log(2*pi)-(n/2)*log(current.res$value)
    current.aic <- 2*(2*ilay-1)-2*current.lik
    current.bic <- log(n)*(2*ilay-1)-2*current.lik
    all.measures[(ilay-1),] <- c(current.err, current.aic, current.bic)
    if(select.measure == "rss"){
      if(current.err < max.err){
        best.res <- current.res
        max.err <- current.err
      }
    }
    else if(select.measure == "aic"){
      if(current.aic < max.aic){
        best.res <- current.res
        max.aic <- current.aic
      }
    }
    else if(select.measure == "bic"){
      if(current.bic < max.bic){
        best.res <- current.res
        max.bic <- current.bic
      }
    }
  }
  best.res$all.measures <- all.measures
  res <- best.res
  return(best.res)
}
#' @title
#' calibrate_joint_nls
#' @description
#' Parameter estimation of more than one VES in a single step.
#' @param ves.list A list of VES objects
#' @param par0 A numeric vector with the initial model
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
#' @family calibration functions
#' @importFrom pracma eye inv
#' @export
calibrate_joint_nls <- function(ves.list, par0, iterations = 30, ireport = 10){
  if(class(ves.list) != "list"){
    stop("ERROR: a list is required as input")
  }
  #
  f_model1 <- function(x, spacing.list){
    nset <- length(spacing.list)
    all.res <- NULL
    npar <- length(x)/nset
    for(iset in 1:nset){
      begin.pos <- 1+(iset-1)*npar
      end.pos <- begin.pos + npar - 1
      par <- 10^x[begin.pos:end.pos]
      spacing <- spacing.list[[iset]]
      res <- apparent_resistivities_simple(par, rves::filt$V1, spacing)
      all.res <- c(all.res, log10(res))
    }
    return(all.res)
  }
  #
  cal.error <- matrix(0.0, nrow = iterations, ncol = 2)
  nves <- length(ves.list)
  if(nves == 1){
    stop('ERROR: At least two VES are required for the joint inversion')
  }
  spacing.list <- list()#ves$ab2
  measured <- NULL
  ndat <- 0
  for(ives in 1:nves){
    spacing.list[[ives]] <- ves.list[[ives]]$ab2
    measured <- c(measured, log10(ves.list[[ives]]$appres))
    ndat <- ndat + length(ves.list[[ives]]$ab2)
  }
  #
  npar <- length(par0)/nves
  nparh <- npar/2
  niter <- iterations
  iter <- 0
  cpar <- log10(par0)
  cpar1 <- as.matrix(cpar, npar, 1)
  mu_max <- 20
  mu <- mu_max
  beta <- 2
  old.par <- NULL
  old.error <- 100
  current.error <- 100
  current.error1 <- 100
  # Initialization of Jacobian and constraint matrices
  # G <- matrix(0.0, nrow = ndat, ncol = length(par0))
  # Cobs <- eye(ndat)
  # CR <- 1*eye((npar*(nves-1)))
  # CCu <- cbind(Cobs, matrix(0.0, nrow = ndat, ncol = (npar*(nves-1)) ))
  # CCl <- cbind(matrix(0.0, nrow = (npar*(nves-1)), ncol = ndat), CR)
  # CC <- rbind(CCu, CCl)
#  R <- matrix(0.0, nrow = (npar*(nves-1)), ncol = length(par0))
#  for(ives in 1:(nves-1)){
#    for(ipar in 1:nparh){
#      begin.pos <- ives*ipar
#      end.pos <- begin.pos + nparh - 1
#      R[end]
#    }
#  }
#  for(ii in 1:(npar*(nves-1))){
#    R[ii,ii] <- 1
#    R[ii,(ii+npar)] <- -1
#  }
  #
  while(iter < niter){
    # Data difference
    cres <- f_model1(cpar, spacing.list)
    d <- as.matrix((cres-measured), length(measured), 1)
    # Calculate jacobian matrix
    J <- jacobian(f_model1, x = cpar, spacing.list = spacing.list)
    #G <- J
    #Gp <- rbind(G, R)
    #delta_r <- -R%*%cpar1
    #dp <- rbind(d, delta_r)
    # calculate parameter correction
    x <- solve(t(J)%*%J+mu*pracma::eye(nves*npar),t(J)%*%d)#/(iter+1)
    #xp <- solve(t(Gp)%*%inv(CC)%*%Gp+mu*eye(nves*npar),t(Gp)%*%inv(CC)%*%dp)
    # calculate error vector
    err <- J%*%x-d
    #errp <- Gp%*%xp-dp
    # Calculate y
    y <- solve(t(J)%*%J+(mu)*eye(nves*npar),t(J)%*%err)
    #yp <- solve(t(Gp)%*%Gp+(mu)*eye(nves*npar),t(Gp)%*%errp)
    E1 <- norm(y)
    #
    old.par <- cpar
    cpar1 <- cpar1 - x #p[1:length(par0)]
    pos.valid <- cpar1 < 0
    if(sum(pos.valid) > 1){
      #cpar1 <- old.par - x
      cpar <- old.par
      print("Correction")
      break
      #print(old.par)
      #print(cpar1)
      #print(xp)
      #print(x)
      #stop('ERROR: the parameter values have become negative')
    }
    cpar <- as.numeric(cpar1)
    cres <- f_model1(cpar, spacing.list)
    pos <- is.nan(cres)
    if(sum(pos) > 1){
      print(cpar)
      stop("ERROR: the forward model evaluated NAN")
    }
    current.error <- mean((cres-measured)**2)
    old.error <- current.error1
    current.error1 <- 100*mean(abs(measured-cres)/measured)
    cal.error[iter+1,1] <- current.error
    cal.error[iter+1,2] <- current.error1
    if(current.error1 > old.error){
      cat("Final step ", iter, ", Error= ", old.error, "\n")
      current.error1 <- old.error
      cpar <- old.par
      break
    }
    if(iter != -1){
      mu2 <- mu/2
      x2 <- solve(t(J)%*%J+mu2*pracma::eye(nves*npar),t(J)%*%d)
      #x2p <- solve(t(Gp)%*%inv(CC)%*%Gp+mu2*pracma::eye(nves*npar),t(Gp)%*%inv(CC)%*%dp)
      err2 <- J%*%x2-d
      #err2p <- Gp%*%x2p-dp
      y2 <- solve(t(J)%*%J+mu2*eye(nves*npar),t(J)%*%err2)
      #y2p <- solve(t(Gp)%*%inv(CC)%*%Gp+mu2*eye(nves*npar),t(Gp)%*%inv(CC)%*%err2p)
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
  #
  # Extract results
  #
  current.res <- list()
  current.rho <- matrix(0.0, nrow = nves, ncol = nparh)
  current.thickness <- matrix(0.0, nrow = nves, ncol = nparh)
  for(ives in 1:nves){
    beginp <- 1 + (ives-1)*npar
    endp <- beginp + npar - 1
    current.par <- 10^cpar[beginp:endp]
    current.res[[ives]] <- current.par
    current.rho[ives,1:nparh] <- current.par[1:nparh]
    current.thickness[ives,1:nparh] <- current.par[(nparh+1):npar]
  }
  #
  res <- list(par = 10^cpar, value = current.error, rel.error = current.error1,
              cal.error = cal.error, rho = current.rho,
              thickness = current.thickness)
  return(res)
}
#' @title
#' calibrate_seq_joint_nls
#' @description
#' Inversion of several VES with the same number of layers using Nonlinear Least-Squares
#' @param ves.list A list of VES objects
#' @param iterations An integer specifying the maximum number of iterations
#' @param ireport An integer specifying the report interval
#' @param max.layers An integer with the maximun number of layers to be tested during the
#' inversion
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
#' @family calibration functions
#' @export
calibrate_seq_joint_nls <- function(ves.list, iterations = 30, ireport = 10,
                                    max.layers = 10){
  if(class(ves.list) != "list"){
    stop('ERROR: A list object is required as input')
  }
  #
  nves <- length(ves.list)
  if(nves == 1){
    stop('ERROR: More than 1 VES is required in the input list')
  }
  current.res <- NULL
  current.par <- NULL
  best.res <- NULL
  current.err <- 0.1
  max.err <- 100
  for(ilay in 2:max.layers){
    def.par <- NULL
    for(ives in 1:nves){
      current.ves <- ves.list[[ives]]
      depth <- current.ves$ab2/2.3
      pos <- round(seq(1, length(depth), length.out = (ilay +1)))
      current.par <- c(rep(mean(current.ves$appres), ilay), diff(depth[pos]))
      def.par <- c(def.par, current.par)
    }
    current.res <- calibrate_joint_nls(ves.list, par0 = def.par,
                                       iterations = iterations,
                                       ireport = ireport)
    current.err <- current.res$rel.error
    if(current.err < max.err){
      best.res <- current.res
      max.err <- current.err
    }
  }
  res <- best.res
  return(best.res)
}
