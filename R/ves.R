#' rves: R package for the interpretation of Vertical Electric Soundings
#'
#' This package includes several functions used in the inversion of the VES.
#'
#' @section base functions:
#' This section includes a description of the main
#'
#' @docType package
#' @name rves
NULL
#' @title
#' ves
#' @description
#' Function to create a VES object
#' @param id a character string defining the VES ID
#' @param ab2 Numeric vector with the values of the electrode spacing
#' @param apprho Numeric vector with the values of the apparent resistivity
#' @return
#' A VES object
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @examples
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' print(sev1a)
#' @export
ves <- function(id = character(0), ab2 = NULL, apprho = NULL){
  if(class(ab2) != 'numeric' | class(apprho) != 'numeric'){
    stop('A numeric vector is required as input')
  }
  #
  ab2a <- NULL
  apprhoa <- NULL
  pos.valid <- ab2 >= 3.0
  if(sum(pos.valid) > 1){
    ab2a <- ab2[pos.valid]
    apprhoa <- apprho[pos.valid]
  }
  else {
    warning("Modifying some of the original data.")
  }
  #
  res <- list(id = id,
              ab2 = ab2a,
              appres = apprhoa,
              interpreted = FALSE,
              rhopar = NULL,
              thickpar = NULL,
              fit.parameters = NULL,
              depth = NULL,
              true.res.direct = NULL,
              true.res.scaling = NULL)
  class(res) <- 'ves'
  invisible(res)
  return(res)
}
#' @title
#' summary.ves
#' @description
#' Function to display a short summary of the apparent resistivity data
#' @param object A VES object
#' @param ... additional parameters for the data.frame summary function
#' @return
#' This function displays on the console a summary of the apparent resistivity data
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @import stringi
#' @export
#' @examples
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' summary(sev1a)
summary.ves <- function(object, ...){
  ves <- object
  cat('Pumping Test: '%s+%ves$id%s+%'\n')
  pdata <- as.data.frame(cbind(ves$ab2, ves$appres))
  names(pdata) <- c("ab/2","Apparent Resistivity(Ohm.m)")
  summary(pdata, ...)
}
#' @title
#' print.ves
#' @description
#' Function to print on the screen a VES object
#' @param x A VES object
#' @param ... Additional parameters to the data.frame print function
#' @return
#' This function prints the information from the VES on the screen
#' @usage \\method{print}{ves}(x, ...)
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @import stringi
#' @export
#' @examples
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' print(sev1a)
print.ves <- function(x, ...){
  ves <- x
  cat('Vertical Electrical Sounding: '%s+%ves$id%s+%'\n')
  pdata <- as.data.frame(cbind(ves$ab2, ves$appres))
  names(pdata) <- c("ab/2","Apparent Resistivity(Ohm.m)")
  print(pdata, ...)
}
#' @title
#' plot.ves
#' @description
#' Function to plot the Vertical Electric Sounding data. This function can create two different
#' types of plots: resistivity and interpretation. The resistivity plot includes the
#' apparent resistivity vs electrode spacing. In the interpretation plot, the values of the
#' apparent resistivity vs electrode spacing and the thickness and real resisitivity of the
#' layers are included.
#' @param x A VES object
#' @param main Title of the plot
#' @param type A character string specifying the plot type. Currently only ves (measurements and
#' earth model), transformation (resistivity-depth transformation) and diagnostic are supported.
#' @param trans.type A character string specifying the transformation type. Only direct, scaling,
#' zohdy and zohdy.smoothed  are currently supported.
#' @param ... Additional parameters to be passed to the plot function
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @family base functions
#' @export
#' @examples
#' library(grid)
#' library(gridExtra)
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' # Conventional VES plot
#' plot(sev1a, type = "ves")
#' # Apply direct transformation to sev1a
#' plot(sev1a, type = "transformation", trans.type = "direct")
#' # Apply scaling transformation to sev1a
#' plot(sev1a, type = "transformation", trans.type = "scaling")
#' # Apply Zohdy transformation to sev1a
#' plot(sev1a, type = "transformation", trans.type = "zohdy")
#' # Apply Zohdy smoothed transformation to sev1a
#' plot(sev1a, type = "transformation", trans.type = "scaling")
#' # Diagnostic plot
#' rho <- c(40,70,30, 20)
#' thick <- c(2,10,50,500)
#' par0 <- c(rho, thick)
#' res.nls <- calibrate_nls(sev1a, par0, iterations = 30, ireport = 5, trace = FALSE)
#' sev1a$rhopar <- res.nls$rho
#' sev1a$thickpar <- res.nls$thickness
#' sev1a$interpreted <- TRUE
#' p1 <- plot(sev1a, type = "ves")
#' print(p1)
plot.ves <- function(x, main = NULL, type = c("ves", "transformation", "diagnostic"),
                     trans.type = c("direct", "scaling", "zohdy", "zohdy.smoothed"),
                     ...){
  if(class(x) != 'ves'){
    stop('A VES object is required as input')
  }
  pres <- NULL
  if(type == "ves"){
    pres <- plot_ves(x, main, ...)
  }
  else if(type == "transformation"){
    pres <- plot_transformation(x, trans.type = trans.type)
  }
  else if(type == "diagnostic"){
    pres <- plot_diagnostic(x)
  }
  return(pres)
}
#' @title
#' plot_ves
#' @description
#' Function to plot the Vertical Electric Sounding data. This function can create two different
#' types of plots: resistivity and interpretation. The resistivity plot includes the
#' apparent resistivity vs electrode spacing. In the interpretation plot, the values of the
#' apparent resistivity vs electrode spacing and the thickness and real resisitivity of the
#' layers are included.
#' @param x A VES object
#' @param main Title of the plot
#' @param ... Additional parameters to be passed to the plot function
#' @importFrom ggplot2 geom_point geom_line scale_x_log10 scale_y_log10 xlab ylab ggtitle
#' @importFrom ggplot2 theme_bw geom_path ggplot aes
#' @importFrom gridExtra grid.arrange
#' @importFrom pracma logseq
#' @export
#' @examples
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' plot_ves(sev1a)
plot_ves <- function(x, main = NULL, ...){
  if(class(x) != 'ves'){
    stop('A VES object is required as input')
  }
  args <- list(...)
  ves <- x
  ves.df <- data.frame(ab2 = ves$ab2, appres = ves$appres)
  pbase <- NULL
  ab2 <- NULL
  appres <- NULL
  res <- NULL
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  yrng <- range(ves$appres)
  log10yrng <- abs(diff(log10(yrng)))
  #print(log10yrng)
  if(!is.null(main)){
    subtitle <- paste(",", main)
  }
  else{
    subtitle <- ""
  }
  if(!ves$interpreted){
    pbase <- ggplot() + geom_point(aes(x = ab2, y = appres), data = ves.df,
                                   col = "red", size = 2) +
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) +
      ylab( expression(paste("Apparent Resitivity ", Omega, phantom() %.% phantom(), "m"))) +
      xlab('AB2(m)') +
      ggtitle(paste("Sounding: ", ves$id, subtitle)) +
      theme_bw()
    if(log10yrng < 1){
      pbase <- pbase + scale_y_log10(breaks = minor_breaks)
    }
    else{
      pbase <- pbase + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)
    }
  }
  if(ves$interpreted){
    if(is.null(ves$rhopar) || is.null(ves$thickpar)){
      stop('ERROR: The interpreted VES has no real resistivities and/or thicknesses.')
    }
    ymn <- 0.8*min(ves$appres, ves$rhopar)
    ymx <- 1.5*max(ves$appres, ves$rhopar)
    #print(c(ymn,ymx))
    fit.parameters <- ves$fit.parameters
    rho <- ves$rhopar
    thick <- ves$thickpar
    spacing <- ves$ab2
    max.depth <- max(spacing/2.3)
    nlayers <- length(rho)
    thick[nlayers] <- max.depth - thick[nlayers-1]
    depth <- 1.0 + cumsum(thick)
    sp <- pracma::logseq(min(spacing), max(spacing), 2*length(ves$ab2))
    res.model <- apparent_resistivities(rho, thick, rves::filt$V1, sp)
    res.model.df <- data.frame(ab2 = res.model$ab2, appres = res.model$appres)
    pos_layers <- vector('numeric', length = (nlayers+1))
    pos_layers[1] <- 1.0
    pos_layers[2:(nlayers+1)] <- depth
    pos_layers1 <- vector('numeric', length = (nlayers*2))
    rho_layers <- vector('numeric', length = (nlayers*2))
    for(ilay in 1:nlayers){
      begin_pos <- (2*ilay-1)
      end_pos <- (2*ilay)
      pos_layers1[begin_pos:end_pos] <- c(pos_layers[ilay], pos_layers[ilay+1])
      rho_layers[begin_pos:end_pos] <- rep(rho[ilay], 2)
    }
    model.par.df <- data.frame(depth = pos_layers1, res = rho_layers)
    pbase <- ggplot() + geom_point(aes(x = ab2, y = appres), data = ves.df,
                                   col = "red", size = 2) +
      geom_line(aes(x = ab2, y = appres), data = res.model.df, col = "blue") +
      geom_path(aes(x = depth, y = res), data = model.par.df, col = "green") +
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) +
      ylab( expression(paste("Apparent Resitivity ", Omega, phantom() %.% phantom(), "m"))) +
      xlab('AB2(m)') +
      ggtitle(paste("Sounding: ", ves$id, subtitle)) +
      theme_bw()
    if(log10yrng < 1){
      pbase <- pbase + scale_y_log10(breaks = minor_breaks)
    }
    else{
      pbase <- pbase + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)
    }
  }
  return(pbase)
}
#' @title
#' plot_transformation
#' @description
#' Function to plot the resistivity-depth transformation from a given VES.
#' @param x A VES object
#' @param trans.type A character string with the type of transformation to apply. Currently only direct and
#' scaling, zohdy and zohdy.smoothed transformation are defined.
#' @return
#' This function returns a ggplot2 object.
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
#' @importFrom ggplot2 geom_point geom_line scale_x_log10 scale_y_log10 xlab ylab ggtitle xlim ylim
#' @examples
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' # Direct Transfor mation
#' plot_transformation(sev1a, trans.type = "direct")
#' # Scaling Transformation
#' plot_transformation(sev1a, trans.type = "scaling")
#' # Zohdy Transformation
#' plot_transformation(sev1a, trans.type = "zohdy")
#' # Zohdy Smoothed Transformation
#' plot_transformation(sev1a, trans.type = "zohdy.smoothed")
plot_transformation <- function(x, trans.type = c("direct", "scaling",
                                                  "zohdy", "zohdy.smoothed")){
  if(class(x) != 'ves'){
    stop('A VES object is required as input')
  }
  # Calculate transformations
  res <- NULL
  ves <- x
  res.trans <- NULL
  if(trans.type == "direct"){
    res.trans <- transform_direct(ves)
  }
  else if(trans.type == "scaling"){
    res.trans <- transform_scaling(ves)
  }
  else if(trans.type == "zohdy"){
    res.trans <- transform_zohdy(ves)
  }
  else if(trans.type == "zohdy.smoothed"){
    res.trans <- transform_smoothed_zohdy(ves)
  }
  else{
    stop("ERROR: Unknown transformation.")
  }
  # Extract results
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  yrng <- range(ves$appres)
  log10yrng <- abs(diff(log10(yrng)))
  depth <- res.trans$depth
  real.res <- res.trans$real.res
  thick1 <- c(diff(depth),500)
  depth1 <- cumsum(thick1)
  depth2 <- log10(res.trans$depth)-0.5*diff(log10(c(0.1,res.trans$depth)))
  res.trans.df <- data.frame(depth = depth, real.res = real.res,
                             depth1 = 10^depth2)
  #
  current.title <- paste0("Resistivity-Depth Transformation(",
                          trans.type,"): Sounding ", ves$id)
  #
  model.par.df <- NULL
  p <- ggplot() + xlab("Depth(m)") +
    ylab( expression(paste("Effective Resitivity ",
                           Omega, phantom() %.% phantom(), "m")) ) +
    ggtitle(current.title)
  ymn <- 0.8*min(ves$appres)
  ymx <- 1.5*max(ves$appres)
  #
  rho <- real.res
  thick <- c(diff(depth),500)
  spacing <- ves$ab2
  nlayers <- length(rho)
  thick[nlayers] <- 1.1*max(depth)
  pos_layers <- vector('numeric', length = (nlayers+1))
  pos_layers[1] <- 0.1
  pos_layers[2:(nlayers+1)] <- depth
  pos_layers1 <- vector('numeric', length = (nlayers*2))
  rho_layers <- vector('numeric', length = (nlayers*2))
  for(ilay in 1:nlayers){
    begin_pos <- (2*ilay-1)
    end_pos <- (2*ilay)
    pos_layers1[begin_pos:end_pos] <- c(pos_layers[ilay], pos_layers[ilay+1])
    rho_layers[begin_pos:end_pos] <- rep(rho[ilay], 2)
  }
  model.par.df <- data.frame(depth = pos_layers1, res = rho_layers)
  ymn <- 0.8*min(ves$appres, ves$rhopar)
  ymx <- 1.5*max(ves$appres, ves$rhopar)
  p <- p + geom_path(aes(x = depth, y = res), data = model.par.df,
                     col = "green")

  #
  x.range <- c(0.1,max(depth))
  y.range <- c(ymn, ymx)
  p <- p + geom_point(aes(x = depth1, y = real.res), data = res.trans.df,
                      color = "blue") +
    theme_bw() +
    scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)
  if(log10yrng < 1){
    p <- p + scale_y_log10(breaks = minor_breaks)
  }
  else{
    p <- p + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)
  }
  #
  return(p)
}
#' @title
#' plot_diagnostic
#' @description
#' Function to create an estimation diagnostic plot composed of:
#' \itemize{
#' \item Scatterplot between the measured and calculated apparement resistivity
#' \item Scatterplot between calculated apparent resistivity and the residuals
#' \item Scatterplot between the calculated apparent resistivity and the standarized residuals
#' \item QQplot of the residuals
#' }
#' @param x A VES object, The VES must be intepreted.
#' @return
#' This function returns a ggplot2 object
#' @author
#' Oscar Garcia-Cabrejo, \email{khaors@gmail.com}
#' @export
#' @importFrom ggplot2 ggplot geom_point geom_smooth coord_equal geom_qq
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @examples
#' data(ves_data1)
#' ab2 <- ves_data1$ab2
#' apprho <- ves_data1$apprho
#' sev1a <- ves(id= "VES1", ab2 = ab2, apprho = apprho)
#' rho <- c(40, 70, 30, 20)
#' thick <- c(2, 10, 50, 500)
#' par0 <- c(rho, thick)
#' res.nls1 <- calibrate_nls(sev1a, par0, iterations = 30, ireport = 5, trace = FALSE)
#' sev1a$rhopar <- res.nls1$rho
#' sev1a$thickpar <- res.nls1$thickness
#' sev1a$interpreted <- TRUE
#' p1 <- plot_diagnostic(sev1a)
#' print(p1)
plot_diagnostic <- function(x){
  if(class(x) != 'ves'){
    stop('A VES object is required as input')
  }
  current.ves <- x
  if(!current.ves$interpreted){
    stop("An interpreted VES is required as input")
  }
  rho <- current.ves$rhopar
  thick <- current.ves$thickpar
  spacing <- current.ves$ab2
  meas.app.rho <- current.ves$appres
  cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                        spacing = spacing)
  residuals.rho <- meas.app.rho-cal.app.rho$appres
  abs.residuals.rho <- sqrt(abs(residuals.rho))
  # Define global variables
  measured <- NULL
  calculated <- NULL
  residuals <- NULL
  abs.residuals <- NULL
  # Plot measured rho vs calculated rho
  df1 <- data.frame(measured = meas.app.rho, calculated = cal.app.rho$appres,
                    residuals = residuals.rho,
                    abs.residuals = abs.residuals.rho)
  p1 <- ggplot() + geom_point(aes(x = measured, y = calculated), data = df1,
                              color = "red") +
    coord_equal() +
    ggtitle("a) Measured vs Calculated")
  #
  p2 <- ggplot(data = df1, aes(x = calculated, y = residuals)) +
    geom_point(color = "red") +
    geom_smooth() +
    ggtitle("b) Residuals")
  #
  p3 <- ggplot(data = df1, aes(x = calculated, y = abs.residuals)) +
    geom_point(color = "red") +
    geom_smooth() +
    ggtitle("c) Absolute Residuals")
  #
  p4 <- ggplot(data = df1, aes(sample = residuals)) + geom_qq(color = "red") +
    ggtitle("d) QQ plot")
  ptot <- arrangeGrob(p1, p2, p3, p4, ncol = 2)
  #ptot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  return(ptot)
}
