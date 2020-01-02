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
              menuItem("Experimental Design", tabName = "ExpDesign", icon = icon('vials'),  startExpanded=TRUE, 
                       menuSubItem("Import design",  tabName = "importExpDesign", selected = TRUE),
                       menuItemOutput("SetUpModel")
                       ),
              menuItemOutput("importData"),
              uiOutput("dataList"),
              menuItemOutput("omics"),
              menuItemOutput("DiffAnalysis"),
              menuItemOutput("CoExpAnalysis")
              
              ),
  tags$br(),
  tags$br(),
  downloadButton("report", "Generate report")
)

body <- dashboardBody(

  tabItems(
    
    #### Import Exp design ####
    ###########################
    

    tabItem(tabName = "importExpDesign",
            fluidRow(
              box(width = 8, status = "warning",

                column(width = 6,
                 # matrix count/abundance input
                 fileInput("Experimental.Design.file", "Import matrix of Experimental Design (txt)",
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                  actionButton("loadExpDesign","load")
                )

              )
            ),
            tags$br(),
            fluidRow(
              uiOutput("ExpDesignTable")
            )
    ),
    
    #### Set Up statistical model & hypothesis ####
    ###############################################
    
    tabItem(tabName = "SetUpModel",
      fluidRow(
        column(width= 12,
          box(status = "warning", width = 12, height = NULL,
              h4("Select the level of reference fo each design factor"),
              fluidRow(
                uiOutput("GetdFactorRef")
              ),
              
              h4("Select the type of the design factor"),
              fluidRow( uiOutput("GetdFactorType") ),
              
              h4("Enter a name for each design factor"),
              fluidRow( uiOutput("GetdFactorName") ),
              
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
    
    #### Import data       ####
    ###########################
    tabItem(tabName = "importData",
            
         box(title = "Load omic data", status = "warning", width = 12, height = NULL,
             h5("[warning] dataset with omics (type == none) was igored..."),
             fluidRow(
                  column(2,
                         # omic type
                         selectInput(inputId='omicType1', label='Omics', 
                                     choices = c("None"="none", "RNAseq"="RNAseq", 
                                                 "Proteomics"="proteomics", "Metabolomics"="Metabolomics"),
                                     selected = "none")
                         ),
                  column(4,
                         
                         # matrix count/abundance input
                         fileInput("data1", "omic count/abundance (Ex.)",
                                   accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                         ),
                  column(4,
                         # metadata/QC bioinfo
                         fileInput("metadataQC1", "QC or metadata (Ex.)",
                                   accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                         ),
                  column(2,
                         # dataset Name
                         textInput(inputId="DataName1", label="Dataset name", value="set1")
                         )
                  
                  #,
                  #column(2, actionButton("removeFactor1", "",
                  #                       icon=icon("times", class = NULL, lib = "font-awesome")))
                  ),
             uiOutput("toAddData2"),
             actionButton("addData","Add data")
             ),
         actionButton("loadData","load")
         ),
    
    #### Data Exploratory & QC ####
    ###############################
    
    tabItem(tabName = "RNAseqExploratoryQC",
            fluidRow(
              box(title = "Raw Data Summary", solidHeader = TRUE, status = "warning", width = 12, height = NULL,
                # library size plot
                column(6, plotOutput("LibSize",   height = "400%")),
                # count distribution plot
                column(6, plotOutput("CountDist", height = "400%"))
              )
            ),
            fluidRow(
              box(title = "Exploratory of Biological and Technical variability", solidHeader = TRUE, width = 12, status = "warning",
                 tabBox( id = "ExplorAnalysisQC", width = 12,

                   tabPanel("Principal component analysis (1/2)",
                        tags$br(),
                        tags$br(),
                        column(width = 2,

                            fluidRow(uiOutput('condColorRaw')),
                            tags$br(),
                            fluidRow(uiOutput('PCA1axisRaw')),
                            fluidRow(uiOutput('PCA2axisRaw')),
                            tags$br(),
                            tags$br(),
                            fluidRow(actionButton("screenshotPCA_QC","Screenshot"))
                            ),
                        column(width = 10,  plotOutput("QCdesignPCARaw"))

                      ),
                   tabPanel("Principal component analysis (2/2)", plotOutput("QCdesignPCA")),
                   tabPanel("Quality check for technical issues", plotOutput("QCdata"))
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

                       #actionButton("RunFiltering","Run Filtering"),
                       verbatimTextOutput("FilterResults")
                  )
                ),
                fluidRow(
                  box( title = "Normalization" , solidHeader = TRUE, status = "warning", width = 12,
                       selectInput(inputId  = "selectNormMethod",
                                   label    = "Method :",
                                   choices  =  list("TMM (edgeR)" = "TMM", "RLE (edgeR)" = "RLE", "upperquartile (edgeR)" = "upperquartile"),
                                   selected = "TMM"),
                       #actionButton("RunNormalization","Run Normalisation")
                  )
                ),
                actionButton("NormValid","Validate")
              ),
              column(8,
                 box(title = "Abundance distribution", solidHeader = TRUE, status = "warning", width = 14 ,  height = NULL,
                      plotOutput("norm.boxplot")
                 )
              )
            ),
            fluidRow(

              box(title = "Principal component analysis", solidHeader = TRUE, status = "warning", width = 12 ,  height = NULL,

                   column(width = 2,
                          fluidRow( uiOutput('condColor') ),
                          tags$br(),
                          fluidRow( uiOutput('PC1axis')),
                          fluidRow( uiOutput('PC2axis')),
                          tags$br(),
                          tags$br(),
                          fluidRow(actionButton("screenshotPCA_Norm","Screenshot"))
                          ),
                   column(width = 10, plotOutput("norm.PCAcoord"))
                  )
              )
            ),
    tabItem(tabName = "ProtExploratoryQC"),
    tabItem(tabName = "ProtProcessing"),
    tabItem(tabName = "MetaExploratoryQC"),
    tabItem(tabName = "MetaProcessing"),
    tabItem(tabName = "DiffAnalysis",
           fluidRow( uiOutput("DiffParam")),
           tags$br(),
           tags$br(),
           fluidRow( uiOutput("ContrastsResults")),
           fluidRow( uiOutput("ResultsMerge"))
    ),
    tabItem(tabName = "CoExpression",
            verbatimTextOutput("Asuivre")
    )
))

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

