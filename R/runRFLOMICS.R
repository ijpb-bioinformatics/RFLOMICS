### ============================================================================
### run interface
### ----------------------------------------------------------------------------
# D. Charif

#' @title Run RFLOMICS interface
#' @description running this function will open the shiny application. 
#' Run the shiny application
#' @param ... More arguments to pass to shinyApp.
#' @return shinyApp
#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinyBS popify
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
# @import shiny 
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' library(RFLOMICS)
#' ## Not run: runRFLOMICS()
#' 
#'
runRFLOMICS <- function(...) {
  # appDir <- system.file("RFLOMICSapp", package = "RFLOMICS")

  # if (appDir == "") {
  #   stop("Could not find example directory. Try re-installing `RFLOMICS`.", call. = FALSE)
  # }

  # shiny::runApp(appDir, display.mode = "normal")
  
  # addResourcePath(prefix = 'docs', system.file('RFLOMICSapp/docs/', package = 'RFLOMICS')) # docs/ n'existe pas ?
  addResourcePath(prefix = 'www', system.file('RFLOMICSapp/www/', package = 'RFLOMICS'))
  shinyApp(ui = rflomicsUI, server = shinyServer(rflomicsServer), ...)
}
