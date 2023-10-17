#' @title runRFLOMICS
#' @description running this function will open the shiny application. 
#' Run the shiny application
#' @return shinyApp
#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinyBS popify
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
# @import shiny 
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom colourpicker colourInput
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' library(RFLOMICS)
#' runRFLOMICS()
#' 
#'
runRFLOMICS <- function(...) {
  # appDir <- system.file("RFLOMICSapp", package = "RFLOMICS")

  # if (appDir == "") {
  #   stop("Could not find example directory. Try re-installing `RFLOMICS`.", call. = FALSE)
  # }

  # shiny::runApp(appDir, display.mode = "normal")
  
  sidebar <- dashboardSidebar(
    
    uiOutput("mysidebar")
  )
  
  body <- dashboardBody({
    
    uiOutput("mycontent")
    
  })
  
  rflomicsUI <- dashboardPage(
    dashboardHeader(title = "RFLOMICS"),
    sidebar,
    body
  ) 
  
  # addResourcePath(prefix = 'docs', system.file('RFLOMICSapp/docs/', package = 'RFLOMICS')) # docs/ n'existe pas ?
  addResourcePath(prefix = 'www', system.file('RFLOMICSapp/www/', package = 'RFLOMICS'))
  shinyApp(ui = rflomicsUI, server = shinyServer(rflomicsServer), ...)
}
