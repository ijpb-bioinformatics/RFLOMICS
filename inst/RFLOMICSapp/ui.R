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
              menuItem("Experimental Design", tabName = "ExpDesign", icon = icon("th"), startExpanded = TRUE,
                       menuSubItem("Import",  tabName = "importExpDesign", selected = TRUE),
                       menuItemOutput("SetUpModel")
                       ),
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
            actionButton("loadExpDesign","load"),
            tags$br(),
            fluidRow(
              tableOutput("ExpDesignTable")
            )
    ),
    tabItem(tabName = "SetUpModel",
      fluidRow(
        column(width= 12, 
          box(status = "warning", width = 12, height = NULL,    
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
            actionButton("ValidF","Valid factor set up")
          )
        )
      ),
      fluidRow(
            column(width= 12, uiOutput("SetModelFormula"))
      ),
      fluidRow(
            column(width= 12, uiOutput("SetContrasts"))
      ),

      tags$br()
    ),
    tabItem(tabName = "importData",
            fluidRow(
                box(title = "RNAseq data", status = "warning", width = 12, height = NULL,
                  column(6,

                     # matrix count/abundance input
                     fileInput("RNAseq.Count.Import.file", "matrix of gene counts as .txt file",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
                  
                  ),
                  column(6,
                     # metadata/QC bioinfo
                     fileInput("RNAseq.QC.Import.file", "QC or metadata as .txt file",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
                  )
                )
            ),
            fluidRow(
              box(title = "Proteomics data",  status = "warning", width = 12, height = NULL,
                column(6,

                     # matrix count/abundance input
                     fileInput("prot.abundances.Import.file", "matrix of protein abundances as .txt file",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                     )
                ),
                column(6,
                     # metadata/QC bioinfo
                     fileInput("prot.QC.Import.file", "QC or metadata as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                              )
                    )
              )
            ),
            fluidRow(
              box(title = "Metabolomics data",  status = "warning", width = 12, height = NULL,
                  
                column(6,

                     # matrix count/abundance input
                     fileInput("metabo.abundances.Import.file", "matrix of metabolite abundances as .txt file",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               )
                     ),
                column(6,
                     # metadata/QC bioinfo
                     fileInput("metabo.QC.Import.file", "QC or metadata as .txt File",
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")
                               )
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
            fluidRow(
              box(title = "Exploratory of Biological and Technical variability", solidHeader = TRUE, width = 12, status = "warning", height = NULL,
                 tabBox( id = "ExplorAnalysisQC", width = 12,
    
                   tabPanel("Quality check of experimental design",  

                        column(width = 2, 
                            fluidRow(uiOutput('condColorRaw')),
                            tags$br(),
                            fluidRow(uiOutput('PCA1axisRaw')),
                            fluidRow(uiOutput('PCA2axisRaw'))
                            #fluidRow(
                            #  radioButtons(inputId  = "PC1raw",
                            #               label    = "Choice of PCs :",
                            #               choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                            #               selected = 1, inline = TRUE),
                            #  
                            #  radioButtons(inputId  = "PC2raw",
                            #               label    = "",
                            #               choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                            #               selected = 2, inline = TRUE)
                            #  )
                        ),
                        column(width = 10,  plotOutput("QCdesignPCARaw"))
                        
                
                      ),
                   tabPanel("Quality check of experimental design",  plotOutput("QCdesignPCA")),
                   tabPanel("Quality check for technical issues",    plotOutput("QCdata"))
                 )
              )
              ),
            tags$br(),
            tags$br()

            ),
    tabItem(tabName = "RNAseqNormalization",
            fluidRow(
              column(4, 
                fluidRow(
                  box( title = "Low Abundance Filtering",  solidHeader = TRUE, status = "warning", width = 12,  
                       numericInput(inputId = "FilterSeuil",
                                    label="Threshold :",
                                    value=0, 0, max=100, 1 ),
                       
                       actionButton("RunFiltering","Run Filtering"),
                       verbatimTextOutput("FilterResults")
                  )
                ),
                fluidRow(
                  box( title = "Normalization" , solidHeader = TRUE, status = "warning", width = 12, 
                       selectInput(inputId  = "selectNormMethod",
                                   label    = "Method :",
                                   choices  =  list("TMM (edgeR)" = "TMM"),
                                   selected = "TMM"),
                       actionButton("RunNormalization","Run Normalisation")
                  )
                )
              ),
              column(8,
                 box(title = "Boxplot", solidHeader = TRUE, status = "warning", width = 14 ,  height = NULL,
                      plotOutput("norm.boxplot")
                 )
              )
            ),
            fluidRow(

              box(title = "PCA", solidHeader = TRUE, status = "warning", width = 12 ,  height = NULL,
                 
                   column(width = 2, 
                          fluidRow( uiOutput('condColor') ),
                          tags$br(),
                          fluidRow( uiOutput('PC1axis')),
                          fluidRow( uiOutput('PC2axis'))
                          ),
                   column(width = 10, plotOutput("norm.PCAcoord"))
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

