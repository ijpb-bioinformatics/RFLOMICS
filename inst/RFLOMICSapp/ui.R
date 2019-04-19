#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinydashboard)

sidebar <- dashboardSidebar(
  sidebarMenu(id="StateSave",
              menuItem("Import Data", tabName = "import",icon = icon('download')),
              menuItemOutput("ExpDesignItem"),
              menuItemOutput("ViewData"),
              menuItemOutput("QualityCheck")
  ),
  downloadButton("report", "Generate report")
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "import",
            fileInput("QC.Import.file", "Import QC data as .txt File",
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")
            ),
            fileInput("Count.Import.file", "Import count matrix as .txt File",
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")
            ),
            actionButton("load","load")
            ),

    tabItem(tabName = "designExp",
            h5("Select the level of reference fo each design factor"),
            fluidRow(
              uiOutput("GetdFactorRef")
              ),
            h5("Select the type of the design factor"),
            fluidRow(
              uiOutput("GetdFactorType")
              ),
            h5("Enter a name for each design factor"),
            fluidRow(
              uiOutput("GetdFactorName")
              ),
            actionButton("ValidF","Valid factor set up"),
            tags$br(),
            tags$br(),
            fluidRow(
            uiOutput("SetModelFormula")
            ),
            fluidRow(
             uiOutput("validM")
            ),
            tags$br(),
            tags$br(),
            fluidRow(
              uiOutput("SetContrasts")
            )
    ),
    tabItem(tabName = "dmatrix", dataTableOutput("colData")
    ),
    tabItem(tabName = "CountMatrix", dataTableOutput("CountMat")
    ),
    tabItem(tabName = "QC",
            numericInput(inputId = "nAxis", label="Number Of PCA axis", value=1,
                         1, max=5, 1)
    ),
    tabItem(tabName = "QCdata",
            plotOutput("QCdata")
    ),
    tabItem(tabName = "QCdesign",
            plotOutput("QCdesign")
    )
  )
)

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

