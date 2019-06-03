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
              menuItemOutput("Filtering"),
              menuItemOutput("ExpDesignItem"),
              menuItemOutput("Exploratory")
  ),
  downloadButton("report", "Generate report")
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "import",
            fluidRow(
              column(6,
                     fileInput("Count.Import.file", "Import matrix of features abundances as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               ),
                     actionButton("load","load")
                     ),
              column(6,
                     fileInput("QC.Import.file", "Import QC data as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               )
                     )
            ),
            tags$br(),
            tags$br(),
            fluidRow(
              column(6,
                    # display report summary
                    box(title = "Data Summary", solidHeader = TRUE, status = "warning", width = 12, tableOutput('SummaryAbundance'))
                    ),
              column(6,
                    box(title = "QC Summary",   solidHeader = TRUE, status = "warning", width = 12, tableOutput('SummaryQC'))
                    )
              ),
            fluidRow(
              column(6,
                    box(title = "Library size", solidHeader = TRUE, status = "warning", width = 12, height = NULL, plotOutput("LibSize", height = "400%"))
                    ),
              column(6,
                    box(title = "Abundance Distribution", solidHeader = TRUE, status = "warning", width = 12, plotOutput("CountDist", height = "400%"))
                    )
              )
            ),
    tabItem(
      tabName = "Filtering",
      h5("Filtering step")
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
    tabItem(tabName = "ExploratoryData", 
            box( status = "warning", width = 12,
                 
                 selectInput("select", label = "Data type",
                             choices = list("Normelized data : TMM" = 1, 
                                            "Unnormalized data" = 2), 
                             selected = 1),
                 tabBox(
                   # The id lets us use input$tabset1 on the server to find the current tab
                   id = "ExplorAnalysis", width = 12,
                   tabPanel("boxplot", plotOutput("norm.boxplot")),
                   tabPanel("QC Design", plotOutput("QCdesign"),
                            numericInput(inputId = "nAxis", label="Number Of PCA axis", value=2, 1, max=5, 1 )
                            ),
                   tabPanel("PCA",
                            fluidRow( plotOutput("norm.PCAcoord")),
                            tags$br(),
                            tags$br(),
                            fluidRow(
                              column(4, uiOutput('condColor')),
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
                            )
                   )
                 )
            ),
    tabItem(tabName = "ExploratoryQC", 
            box( status = "warning", width = 12,
                 tabPanel("QC Data",   plotOutput("QCdata"))
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

