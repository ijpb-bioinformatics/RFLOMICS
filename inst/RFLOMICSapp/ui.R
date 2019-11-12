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
              menuItem("Import ExpDesign", tabName = "importExpDesign",icon = icon('download')),
              menuItemOutput("ExpDesignItem"),
              menuItemOutput("importData"),
              menuItemOutput("omics"),
              #menuItemOutput("RNAseq") ,
              #menuItemOutput("metabolome") ,
              #menuItemOutput("proteome"),
              menuItemOutput("Exploratory"),
              menuItemOutput("DiffAnalysis")
              #menuItemOutput("CoExpAnalysis")
              ),
  tags$br(),
  tags$br(),
  downloadButton("report", "Generate report")
)

body <- dashboardBody(

  tabItems(
    tabItem(tabName = "importExpDesign",
            fluidRow(
              column(6,

                     # matrix count/abundance input
                     fileInput("Experimental.Design.file", "Import matrix of Experimental Design (txt)",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
              )
            ),
            actionButton("loadExpDesign","load")
    ),
    tabItem(tabName = "importData",
            fluidRow(
              column(6,

                     # matrix count/abundance input
                     fileInput("RNAseq.Count.Import.file", "RNAseq : Import matrix of gene counts (txt)",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
              ),
              column(6,
                     # metadata/QC bioinfo
                     fileInput("RNAseq.QC.Import.file", "Import RNAseq QC/metadata (txt)",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
              )
            ),
            fluidRow(
              column(6,

                     # matrix count/abundance input
                     fileInput("prot.abundances.Import.file", "Proteom : Import matrix of protein abundances (txt)",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
              ),
              column(6,
                     # metadata/QC bioinfo
                     fileInput("prot.QC.Import.file", "Import prot QC/metadata as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                              )
                    )
            ),
            fluidRow(
              column(6,

                     # matrix count/abundance input
                     fileInput("metabo.abundances.Import.file", "Metabo : Import matrix of metabolite abundances (txt)",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               )
                     ),
              column(6,
                     # metadata/QC bioinfo
                     fileInput("metabo.QC.Import.file", "Import meta QC/metadata as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               )
                     )
              ),
            actionButton("loadData","load")#,

            #column(6,
            #       # display metadata count summary
            #       box(title = "QC/metadata Summary",   solidHeader = TRUE, status = "warning", width = 12,
            #           tableOutput('SummaryQC'))
            #)

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
             uiOutput("validModelFormula")
            ),
            tags$br(),
            tags$br(),
            fluidRow(
              uiOutput("SetContrasts")
            ),
            fluidRow(
              uiOutput("validContrasts")
            ),
            tags$br()
    ),
    tabItem(tabName = "RNAseqExploratoryQC",
            fluidRow(
              box(title = "Raw Data Summary", solidHeader = TRUE, status = "warning", width = 12, height = NULL,
                column(6,
                    # library size plot
                    plotOutput("LibSize", height = "400%")
                ),
                column(6,
                    # count distribution plot
                     plotOutput("CountDist", height = "400%")
                )
              )
            ),
            box(title = "Exploratory of Biological and Technical variability", solidHeader = TRUE, status = "warning", width = 12,height = NULL,
                 tabBox(
                   id = "ExplorAnalysisQC", width = 12,
                   tabPanel("Quality check of experimental design",  plotOutput("QCdesign")),
                   tabPanel("Quality check for technical issues",    plotOutput("QCdata"))
                 )
            ),
            tags$br(),
            tags$br()
            #box( status = "warning", width = 12,
                 #tabBox(
                   # The id lets us use input$tabset1 on the server to find the current tab
                   #id = "ExplorAnalysis", width = 12,
                   # summay table of filtering step
                   #tabPanel("summary",
                   #
                  #          column(6,
                  #
                  #          ),
                  #          column(6,
                  #                 # display matrix count summary
                  #
                  #          ),
                  #          column(6,
                  #                 # display matrix count summary
                  #
                  #          )
                  # )#,
                  # tabPanel("PCA",
                  #          selectInput("selectData2", label = "Data :",
                  #                      choices = list("Normelized data" = "norm",
                  #                                     "Unnormalized data" = "raw"),
                  #                      selected = "norm"),

                  #          fluidRow( plotOutput("norm.PCAcoord") ),
                  #          tags$br(),
                  #          tags$br(),
                  #          fluidRow(
                  #            column(width = 6, uiOutput('condColor')),
                  #            column(width = 6,
                  #                   radioButtons(inputId  = "PC1",
                  #                                label    = "Choice of PCs :",
                  #                                choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                  #                                selected = 1, inline = TRUE),
                  #
                  #                   radioButtons(inputId  = "PC2",
                  #                                label    = "",
                  #                                choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                  #                                selected = 2, inline = TRUE))
                  #
                  #            )
                  #          )
                  # )
                 #)

            ),
    tabItem(tabName = "RNAseqNormalization",
            fluidRow(
              box( title = "Low Abundance Filtering",  solidHeader = TRUE, status = "warning", width = 6,
                   numericInput(inputId = "FilterSeuil",
                                label="Threshold :",
                                value=0, 0, max=100, 1 ),
                   actionButton("RunFiltering","Run Filtering")
              ),
              box( title = "Normalization" , solidHeader = TRUE, status = "warning", width = 6,
                   selectInput(inputId  = "selectNormMethod",
                               label    = "Method :",
                               choices  =  list("TMM (edgeR)" = "TMM"),
                               selected = "TMM"),
                   actionButton("RunNormalization","Run Normalisation")
              )
            ),
            fluidRow(

              box(title = "Processed Data Summary", solidHeader = TRUE, status = "warning", width = 12 ,
                  column(4,
                         tableOutput('SummaryAbundance')
                  ),
                  column(8,
                         plotOutput("norm.boxplot")
                  )
              )
            )
    ),
   tabItem(tabName = "DiffAnalysis",
           fluidRow(
             column(5,
                           selectInput("AnaDiffMethod", label = "Method :",
                                       choices = list("glmfit (edgeR)"="edgeRglmfit",
                                                      "limma" = "limma"),
                                       selected = "edgeRglmfit")
           ),
            column(3,
                            numericInput(inputId = "FDRSeuil",
                                    label="FDR :",
                                    value=0.05, 0, max=1, 0.01)
            ),
           column(5,
                  actionButton("runAnaDiff","Run the differential analysis")
           )),
           tags$br(),
           tags$br(),
           fluidRow(
             column(width = 8,
             uiOutput("ContrastsResults")
             )
           )
    #tabItem(tabName = "CoExpAnalysis")
    )
))

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

