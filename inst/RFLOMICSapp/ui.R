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
              menuItemOutput("QualityCheck"),
              menuItemOutput("Normalization")
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
            fileInput("Count.Import.file", "Import matrix of features abundances as .txt File",
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
    tabItem(tabName = "AbundanceMatrix", dataTableOutput("CountMat")
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
    ),
    tabItem(tabName = "Norm", 
            
            box(title = "Normalization", solidHeader = TRUE, status = "warning", width = 12,
            
                tabBox(
                  # The id lets us use input$tabset1 on the server to find the current tab
                  id = "tabset1", width = 12,
                  tabPanel("boxplot", plotOutput("norm.boxplot")),
                  tabPanel("PCA 1/2", 
                           fluidRow(
                             plotOutput("norm.PCAbar")
                           ),
                           tags$br(),
                           tags$br(),
                           fluidRow(
                             column(width = 4,
                                    uiOutput('condColor1')
                             )
                           )
                           
                           ),
                  tabPanel("PCA 2/2",     

                           fluidRow(
                             plotOutput("norm.PCAcoord")
                             ),
                           tags$br(),
                           tags$br(),
                           fluidRow(
                             column(width = 4,
                                    uiOutput('condColor')
                                    ),
                             column(width = 4,
                                    radioButtons(inputId = "PC1", 
                                                 label = "Choice of PCs :", 
                                                 choices = list("PC1" = 1, 
                                                                "PC2" = 2, 
                                                                "PC3" = 3),
                                                 selected = 1, inline = TRUE),
                                    radioButtons(inputId = "PC2", 
                                                 label = "", 
                                                 choices = list("PC1" = 1, 
                                                                "PC2" = 2, 
                                                                "PC3" = 3),
                                                 selected = 2, inline = TRUE)
                                    )
                             )
                           ),
                  tabPanel("MA-plot", textOutput("NormText"))
                )
            )        
    )
  )
)

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

