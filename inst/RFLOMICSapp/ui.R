#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinyFiles)

sidebar <- dashboardSidebar(

  sidebarMenu(id="StateSave",
              menuItem(text = "Welcome", tabName = "coverPage", icon = icon('dna')),
              menuItem(text = "Load Data", tabName = "importData", icon = icon('download')),
              menuItemOutput(outputId = "SetUpModelMenu"),

              # menuItem(text = "Experimental Design", tabName = "ExpDesign", icon = icon('vials'),
              #          menuSubItem(text = "Import design",  tabName = "importExpDesign"),
              #          menuItemOutput(outputId = "SetUpModelMenu")
              #          ),
              menuItemOutput(outputId = "omics"),
              menuItemOutput(outputId = "Integration")

              ),
  tags$br(),
  tags$br(),
  uiOutput("runReport")
  #downloadButton(outputId = "report", label = "Generate report")
)

body <- dashboardBody({

  items <- c(

    list(

        #### Cover Page        ####
        ###########################

        tabItem(tabName = "coverPage",
                coverPage

               # shinyDirButton(id = "dir0", label = "Input directory", title = "Upload"),
        ),

        #### Import data       ####
        ###########################
        tabItem(tabName = "importData",

                LoadOmicsDataUI("data")
        ),

        #### Set Up statistical model & hypothesis ####
        ###############################################
        tabItem(tabName = "SetUpModel",

                GLM_modelUI("model")
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

                       QCNormalizationTabUI(paste0("RNAseq",i))

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

                       QCNormalizationTabUI(paste0("proteomics",j))
              ),
              #### Diff analysis  ####
              ######################################
              tabPanel("Diff Prot Expression",
                       tags$br(),
                       tags$br(),
                       DiffExpAnalysisUI(paste0("proteomics",j))
              ),
              #### Co-expression analysis  ####
              ######################################
              tabPanel("Proteins CoExpression",
                       tags$br(),
                       tags$br(),
                       CoSeqAnalysisUI(paste0("proteomics",j))
                       #verbatimTextOutput("Asuivre")
              ),
              #### enrichment analysis  ####
              ######################################
              tabPanel("Annotation Enrichment",
                       tags$br(),
                       tags$br(),
                       AnnotationEnrichmentUI(paste0("proteomics",j))
              )
            )
          )
    }),
  #### metabolome Analysis  #######
  ###############################
  lapply(1:10, function(k){

    tabItem(tabName = paste0("MetaAnalysis", k),
            tabsetPanel(

              #### Data Exploratory & QC ####
              ###############################
              tabPanel("Data Exploratory",
                       tags$br(),
                       tags$br(),

                       QCNormalizationTabUI(paste0("metabolomics",k))
              ),#### Diff analysis  ####
              ######################################
              tabPanel("Diff Metabo Expression",
                       tags$br(),
                       tags$br(),
                       DiffExpAnalysisUI(paste0("metabolomics",k))
              ),
              #### Co-expression analysis  ####
              ######################################
              tabPanel("Metabolites CoExpression",
                       tags$br(),
                       tags$br(),
                       CoSeqAnalysisUI(paste0("metabolomics",k))
                       #verbatimTextOutput("Asuivre")
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

