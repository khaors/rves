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
#' @export
ves <- function(id = character(0), ab2 = NULL, apprho = NULL){
  if(class(ab2) != 'numeric' | class(apprho) != 'numeric'){
    stop('A numeric vector is required as input')
  }
  res <- list(id = id, ab2 = ab2, appres = apprho,
              interpreted = FALSE,
              rhopar = NULL,
              thickpar = NULL,
              fit.parameters = NULL)
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
#' @param ... Additional parameters to be passed to the plot function
#' @importFrom ggplot2 geom_point geom_line scale_x_log10 scale_y_log10 xlab ylab ggtitle
#' @importFrom ggplot2 theme_bw geom_path ggplot aes
#' @importFrom gridExtra grid.arrange
#' @importFrom pracma logseq
#' @export
plot.ves <- function(x, ...){
  if(class(x) != 'ves'){
    stop('A VES object is required as input')
  }
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
  if(!ves$interpreted){
    pbase <- ggplot() + geom_point(aes(x = ab2, y = appres), data = ves.df,
                                   col = "red", size = 2) +
      scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) +
      ylab( expression(paste("Apparent Resitivity ", Omega, phantom() %.% phantom(), "m"))) +
      xlab('AB2(m)') +
      ggtitle(paste("Profile: ", ves$id)) +
      theme_bw()
    if(log10yrng < 1){
      pbase <- pbase + scale_y_log10(breaks = minor_breaks)
    }
    else{
      pbase <- pbase + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)
    }
    print(pbase)
  }
  if(ves$interpreted){
    ymn <- 0.8*min(ves$appres, ves$rhopar)
    ymx <- 1.5*max(ves$appres, ves$rhopar)
    #print(c(ymn,ymx))
    fit.parameters <- ves$fit.parameters
    rho <- ves$rhopar
    thick <- ves$thickpar
    spacing <- ves$ab2
    nlayers <- length(rho)
    depth <- min(spacing) + cumsum(thick)
    sp <- pracma::logseq(min(spacing), max(spacing), 2*length(ves$ab2))
    res.model <- apparent_resistivities(rho, thick, rves::filt$V1, sp)
    res.model.df <- data.frame(ab2 = res.model$ab2, appres = res.model$appres)
    pos_layers <- vector('numeric', length = (nlayers+1))
    pos_layers[1] <- min(spacing)
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
      ggtitle(paste("Survey: ", ves$id)) +
      theme_bw()
    if(log10yrng < 1){
      pbase <- pbase + scale_y_log10(breaks = minor_breaks)
    }
    else{
      pbase <- pbase + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)
    }
    print(pbase)
  }
  return(pbase)
}

