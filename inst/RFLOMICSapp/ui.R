#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinydashboard)
# source("DataExploratoryModules.R")
# source("NormalizationModules.R")
# source("DiffExpressionModules.R")
# source("CoExpressionModules.R")
# source("commonModules.R")

sidebar <- dashboardSidebar(
  sidebarMenu(id="StateSave",
              menuItem("Experimental Design", tabName = "ExpDesign", icon = icon('vials'),  startExpanded=TRUE, 
                       menuSubItem("Import design",  tabName = "importExpDesign", selected = TRUE),
                       menuItemOutput("SetUpModel")
                       ),
              menuItemOutput("importData"),
              #uiOutput("dataList"),
              menuItemOutput("omics"),
              menuItemOutput("DiffAnalysis"),
              menuItemOutput("CoExpAnalysis")
              
              ),
  tags$br(),
  tags$br(),
  downloadButton("report", "Generate report")
)

body <- dashboardBody({

  items <- c( 
    
    list(
  
  #tabItems(
    
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
         )),
    
    #### RNAseq Analysis  ####
    ###############################
    lapply(1:10, function(i){ 
      tabItem(tabName = paste0("RNAseqAnalysis", i),
            tabsetPanel(
              
              #### Data Exploratory & QC ####
              ###############################
              tabPanel("Data Exploratory",
                       tags$br(),
                       tags$br(),
                       
                       RNAseqDataExplorTabUI(i)
                       
                       ),
              
              #### Data Filter & Normalisation  ####
              ######################################
              tabPanel("Filter & Normalization",
                       tags$br(),
                       tags$br(),
                       RNAseqDataNormTabUI(i)
                       ),
              #### Diff analysis  ####
              ######################################
              tabPanel("Diff Gene Expression",
                       tags$br(),
                       tags$br(),
                       DiffExpAnalysisUI(i)
                       ),
              #### Co-expression analysis  ####
              ######################################
              tabPanel("Gene CoExpression", 
                       #verbatimTextOutput("Asuivre")
              ),
              tabPanel("Annotation Enrichment", 
                       #verbatimTextOutput("Asuivre")
              )
              
            )
    )
  })
    #### RNAseq Analysis  ####
    ###############################
    #tabItem(tabName = "proteomicsAnalysisSteps",
    #        tabsetPanel(
    #          
    #          #### Data Exploratory & QC ####
    #          ###############################
    #          
    #          tabPanel("ProtExploratoryQC"),
    #          tabPanel("ProtProcessing"),
    #          tabPanel("CoExpression"),
    #          tabPanel("Annotation")
    #        )
    #)
 )
do.call(tabItems, items)
})

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

