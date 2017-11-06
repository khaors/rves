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
#' @importFrom grDevices colors
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @export
plot.ves <- function(x, ...){
  if(class(x) != 'ves'){
    stop('A VES object is required as input')
  }
  ves <- x
  if(!ves$interpreted){
    plot(ves$ab2, ves$appres, type = 'p', log = 'xy', xlab = 'ab2(m)',
         ylab = 'Apparent Resistivity(Ohm.m)', main = ves$id)
  }
  if(ves$interpreted){
    ymn <- 0.8*min(ves$appres, ves$rhopar)
    ymx <- 1.5*max(ves$appres, ves$rhopar)
    plot(ves$ab2, ves$appres, type = 'p', log = 'xy', xlab = 'ab2(m)',
         ylab = 'Apparent Resistivity(Ohm.m)', ylim = c(ymn, ymx), main = ves$id)
    fit.parameters <- ves$fit.parameters
    rho <- ves$rhopar
    thick <- ves$thickpar
    nlayers <- length(rho)
    depth <- cumsum(thick)
    pos_layers <- vector('numeric', length = (nlayers+1))
    pos_layers[1] <- 1
    pos_layers[2:(nlayers+1)] <- depth
    pos_layers1 <- vector('numeric', length = (nlayers*2))
    rho_layers <- vector('numeric', length = (nlayers*2))
    for(ilay in 1:nlayers){
      begin_pos <- (2*ilay-1)
      end_pos <- (2*ilay)
      pos_layers1[begin_pos:end_pos] <- c(pos_layers[ilay], pos_layers[ilay+1])
      rho_layers[begin_pos:end_pos] <- rep(rho[ilay], 2)
    }
    lines(pos_layers1, rho_layers, col = 'red')
    res.model <- apparent_resistivities(rho,thick, as.matrix(rves::filt$V1))
    lines(res.model$ab2, res.model$appres, col = "blue")
  }
}

