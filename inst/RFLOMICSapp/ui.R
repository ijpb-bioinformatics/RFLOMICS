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
              menuItemOutput("Exploratory") ,
              menuItemOutput("DiffAnalysis"),
              menuItemOutput("CoExpAnalysis")
              ),
  downloadButton("report", "Generate report")
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "import",
            fluidRow(
              column(6,
                     # matrix count/abundance input
                     fileInput("Count.Import.file", "Import matrix of features abundances as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               ),
                     actionButton("load","load")
                     ),
              column(6,
                     # metadata/QC bioinfo
                     fileInput("QC.Import.file", "Import QC/metadata as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               ),
                     selectInput("ExperimentType", label = "Data :",
                                                 choices = list("RNA-seq" = "RNAseq"), 
                                                selected = "RNAseq")
                     )
            ),
            tags$br(),
            tags$br(),
            fluidRow(
              column(6,
                    # display matrix count summary
                    box(title = "Data Summary", solidHeader = TRUE, status = "warning", width = 12, 
                        tableOutput('SummaryAbundance'))
                    ),
              column(6,
                    # display metadata count summary
                    box(title = "QC/metadata Summary",   solidHeader = TRUE, status = "warning", width = 12, 
                        tableOutput('SummaryQC'))
                    )
              ),
            fluidRow(
              column(6,
                    # library size plot
                    box(title = "Library size", solidHeader = TRUE, status = "warning", width = 12, height = NULL, 
                        plotOutput("LibSize", height = "400%"))
                    ),
              column(6,
                     # count distribution plot
                    box(title = "Abundance Distribution", solidHeader = TRUE, status = "warning", width = 12, 
                        plotOutput("CountDist", height = "400%"))
                    )
              )
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
            box( title = "Filtering" ,      solidHeader = TRUE, status = "warning", width = 6,
                 numericInput(inputId = "FilterSeuil", 
                              label="Threshold :", 
                              value=0, 0, max=10, 1 )),
            
            box( title = "Normalization" , solidHeader = TRUE, status = "warning", width = 6, 
                 selectInput(inputId  = "selectNormMethod", 
                             label    = "Method :",
                             choices  =  list("TMM (edgeR)" = "TMM"), 
                             selected = "TMM")),
            actionButton("Norm","Update"),
            tags$br(),
            tags$br(),
            box( status = "warning", width = 12,
                 tabBox(
                   # The id lets us use input$tabset1 on the server to find the current tab
                   id = "ExplorAnalysis", width = 12,
                   # summay table of filtering step
                   tabPanel("summary",
                            column(width=4,
                                   tags$br(),
                                   tags$br(),
                                   tableOutput('SummaryAbundance2'),
                                   textOutput("SummaryText")
                                   ),
                            column(width=8,
                                   tags$br(),
                                   tags$br(),
                                   plotOutput("NormFact")
                                   )
                            ),
                   tabPanel("boxplot", plotOutput("norm.boxplot")),
                   tabPanel("PCA",
                            selectInput("selectData2", label = "Data :",
                                        choices = list("Normelized data" = "norm", 
                                                       "Unnormalized data" = "raw"), 
                                        selected = "norm"),

                            fluidRow( plotOutput("norm.PCAcoord") ),
                            tags$br(),
                            tags$br(),
                            fluidRow(
                              column(width = 6, uiOutput('condColor')),
                              column(width = 6,
                                     radioButtons(inputId  = "PC1", 
                                                  label    = "Choice of PCs :", 
                                                  choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                                                  selected = 1, inline = TRUE),
                                  
                                     radioButtons(inputId  = "PC2", 
                                                  label    = "", 
                                                  choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                                                  selected = 2, inline = TRUE))
                              
                              )
                            )
                   )
                 )
            ),
    tabItem(tabName = "ExploratoryQC", 
            box( status = "warning", width = 12,
                 selectInput("selectData", label = "Data :",
                                           choices = list("Normalized data" = "norm", 
                                                          "Unnormalized data" = "raw"), 
                                           selected = "norm"),
                 tabBox(
                   id = "ExplorAnalysisQC", width = 12,
                   tabPanel("QC Design",  plotOutput("QCdesign")),
                   tabPanel("QC Data",    plotOutput("QCdata"))
                 )
            )
    ),
    tabItem(tabName = "DiffAnalysis"),
    tabItem(tabName = "CoExpAnalysis")
    )
)

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

