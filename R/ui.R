#' @title rflomicsUI
#' @export
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
#' @rawNamespace import(shiny, except = renderDataTable)
rflomicsUI <- function(){
   uiout <- dashboardPage(
      dashboardHeader(title = "RFLOMICS"),
      dashboardSidebar(
        
        uiOutput("mysidebar")
      ),
      dashboardBody({
        
        uiOutput("mycontent")
        
      })
    )
   return(uiout)
  }

rflomicsUI()

