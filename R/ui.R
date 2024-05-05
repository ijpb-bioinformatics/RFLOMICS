#' @title rflomicsUI
#' @keywords internal
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
#' @rawNamespace import(shiny, except = renderDataTable)
#' @noRd
#' @return a user interface.
rflomicsUI <- function(){
  
  ui <- tagList(
    
    tags$head(
      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = "www/style.css")
    ),
    
    dashboardPage(
      
      dashboardHeader(title = "RFLOMICS"),
      dashboardSidebar(
        
        uiOutput("mysidebar")
      ),
      dashboardBody({
        
        uiOutput("mycontent")
      })
    ))
  return(ui)
}

# rflomicsUI()

