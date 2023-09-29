#' @title runRFLOMICS
#' @description running this function will open the shiny application. 
#' Run the shiny application
#' @return shinyApp
#' @export
#'
runRFLOMICS <- function() {
  appDir <- system.file("RFLOMICSapp", package = "RFLOMICS")

  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `RFLOMICS`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
