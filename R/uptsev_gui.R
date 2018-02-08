#' @title
#' uptsev_gui
#' @description
#' uptsev GUI
#' @import shiny
#' @author
#' Oscar Garcia-Cabrejo \email{khaors@gmail.com}
#' @export
uptsev_gui <- function() {
  appDir <- system.file("Shiny", "uptsev", package = "rves")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `rves`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

