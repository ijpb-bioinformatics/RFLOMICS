#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinydashboard)
#library(shinyjs)
# source("DataExploratoryModules.R")
# source("NormalizationModules.R")
# source("DiffExpressionModules.R")
# source("CoExpressionModules.R")
# source("commonModules.R")

sidebar <- dashboardSidebar(
  
  sidebarMenu(id="StateSave",
              menuItem("Presentation", tabName = "coverPage", icon = icon('dna'), startExpanded=TRUE, selected = TRUE),
              menuItem("Experimental Design", tabName = "ExpDesign", icon = icon('vials'), 
                       menuSubItem("Import design",  tabName = "importExpDesign"),
                       menuItemOutput("SetUpModelMenu")
                       ),
              menuItemOutput("importData"),
              menuItemOutput("omics"),
              menuItemOutput("Integration")
              
              ),
  tags$br(),
  tags$br(),
  downloadButton("report", "Generate report")
)

body <- dashboardBody({

  items <- c( 
    
    list(
  
        #### Cover Page        ####
        ###########################
        
        tabItem(tabName = "coverPage",
                verbatimTextOutput("print")
              ),
      
        #### Import Exp design ####
        ###########################
        
    
        tabItem(tabName = "importExpDesign",
                
              ExperimentalDesignUI("Exp"),
        ),
        
        #### Set Up statistical model & hypothesis ####
        ###############################################
        
        tabItem(tabName = "SetUpModel",
              
                GLM_modelUI("model")
        ),
        
        #### Import data       ####
        ###########################
        tabItem(tabName = "importData",
                
             LoadOmicsDataUI("data")
        )
    ),
    
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
                       
                       RNAseqDataExplorTabUI(paste0("RNAseq",i))
                       
                       ),
              
              #### Data Filter & Normalisation  ####
              ######################################
              tabPanel("Filter & Normalization",
                       tags$br(),
                       tags$br(),
                       RNAseqDataNormTabUI(paste0("RNAseq",i))
                       ),
              #### Diff analysis  ####
              ######################################
              tabPanel("Diff Gene Expression",
                       tags$br(),
                       tags$br(),
                       DiffExpAnalysisUI(paste0("RNAseq",i))
                       ),
              #### Co-expression analysis  ####
              ######################################
              tabPanel("Gene CoExpression", 
                       tags$br(),
                       tags$br(),
                       CoSeqAnalysisUI(paste0("RNAseq",i))
                       #verbatimTextOutput("Asuivre")
                       ),
              #### enrichment analysis  ####
              ######################################
              tabPanel("Annotation Enrichment",
                       tags$br(),
                       tags$br(),
                       AnnotationEnrichmentUI(paste0("RNAseq",i))
              )
              
            )
    )
      
  }),
  
  #### Proteome Analysis  #######
  ###############################
  lapply(1:10, function(j){ 
    
    tabItem(tabName = paste0("ProtAnalysis", j),
            tabsetPanel(
              
              #### Data Exploratory & QC ####
              ###############################
              tabPanel("Data Exploratory",
                       tags$br(),
                       tags$br(),
                       
                       ProtMetaDataExplorTabUI(paste0("proteomics",j))
              )
            )
          )
    }),
  #### metabolome Analysis  #######
  ###############################
  lapply(1:10, function(i){ 
    
    tabItem(tabName = paste0("MetaAnalysis", i),
            tabsetPanel(
              
              #### Data Exploratory & QC ####
              ###############################
              tabPanel("Data Exploratory",
                       tags$br(),
                       tags$br(),
                       
                       ProtMetaDataExplorTabUI(paste0("metabolomics",i))
              )
            )
    )
  })
 )
do.call(tabItems, items)
})

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "RFLOMICS"),
  sidebar,
  body
)

