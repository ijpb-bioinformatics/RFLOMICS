#' @title rflomicsUI
#' @export
#' @importFrom shinydashboard dashboardSidebar dashboardBody dashboardHeader dashboardPage
#' @rawNamespace import(shiny, except = renderDataTable)
rflomicsUI <- dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  dashboardSidebar(
    
    uiOutput("mysidebar")
  ),
  dashboardBody({
    
    uiOutput("mycontent")
    
  })
) 

