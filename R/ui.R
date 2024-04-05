#' @title rflomicsUI
#' @keywords internal
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
#' @rawNamespace import(shiny, except = renderDataTable)
#' @noRd
#' @return a user interface.
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

# rflomicsUI()

