#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinyBS popify
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput 
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom colourpicker colourInput
#' @importFrom magrittr "%>%"
#' @importFrom rmarkdown render
rflomicsServer <- function(input, output, session) {
  
  #### Increasing maximum possible size of loaded files (default is only 5MB)
  # https://stackoverflow.com/questions/18037737/how-to-change-maximum-upload-size-exceeded-restriction-in-shiny-and-save-user
  options(shiny.maxRequestSize = 100*1024^2) # 100 MB limit.
  
  # This is to get the desired menuItem selected initially.
  # selected=T seems not to work with a dynamic sidebarMenu.
  observeEvent(session, {
    updateTabItems(session = session, inputId = "tabs", selected = "coverPage")
  })
  
  #############################################
  # reactive value for reinitialisation of UIoutput
  rea.values <- reactiveValues(
    validate.status = 0,
    loadData = FALSE,
    model    = FALSE,
    analysis = FALSE,
    resetAna = FALSE,
    report   = FALSE,
    
    exampleData = NULL,
    
    datasetList   = NULL,
    contrastList  = NULL,
    Contrasts.Sel = NULL,
    datasetProcess  = NULL,
    datasetDiff     = NULL, # list of dataset names with diff results
    datasetCoEx     = NULL,
    datasetDiffAnnot= NULL,
    datasetCoExAnnot= NULL
  )
  
  #############################################
  # dynamic sidebar menu #
  output$mysidebar <- renderUI({
    
    tagList(
      sidebarMenu(id="tabs",
                  menuItem(text = "Welcome", tabName = "coverPage", icon = icon('dna'), selected = TRUE),
                  menuItem(text = "Glossary page", tabName = "GlossaryPage", icon = icon("address-book")),
                  menuItem(text = "Load Data", tabName = "importData", icon = icon('download')),
                  menuItemOutput(outputId = "SetUpModelMenu"),
                  menuItemOutput(outputId = "omics"),
                  menuItemOutput(outputId = "Integration")
      ),
      
      tags$br(),
      tags$br(),
      uiOutput("runReport"),
      tags$br(),
      tags$br(),
      uiOutput("downloadResults")
      
    )
  })
  
  
  #### Item for each omics #####
  # display omics Item
  output$omics <- renderMenu({
    
    validate({
      need(rea.values$analysis == TRUE, message="")
    })
    
    menuItem(text = "Omics Analysis", tabName = "OmicsAnalysis", icon = icon('chart-area'),
             list(
               lapply(names(rea.values$datasetList), function(omics){
                 
                 lapply(names(rea.values$datasetList[[omics]]), function(i){
                   menuSubItem(text = rea.values$datasetList[[omics]][[i]],
                               tabName = paste0(omics, "Analysis", i))
                 })
               }),
               menuItemOutput("omicsSumUI")
             )
    )
  })
  
  output$omicsSumUI <- renderMenu({
    
    validate(
      need(rea.values$analysis == TRUE && length(rea.values$datasetProcess) >= 2, message = "")
    )
    
    menuSubItem(text = "Summary of Analyses", tabName = "omicsSum" )
  })
  
  #### Item for each data integration tools #####
  # display tool Item
  output$Integration <- renderMenu({
    
    validate({
      need(rea.values$analysis == TRUE && length(rea.values$datasetProcess) >= 2, message = "")
    })
    
    menuItem(text = "Data Integration", tabName = "OmicsIntegration", icon = icon('network-wired'), startExpanded = FALSE,selected = FALSE,
             menuSubItem(text = "with MOFA", tabName = "withMOFA" ),
             menuSubItem(text = "with MixOmics", tabName = "withMixOmics")
    )
  })
  
  #### Item for report #####
  output$runReport <- renderUI({
    if(is.null(rea.values$datasetProcess) || 
       length(rea.values$datasetProcess) != length(rea.values$datasetList)) return()

    column(
      width = 12, 
      downloadButton(outputId = "report", 
                     label = "Generate report", class = "butt"),
      tags$head(
        tags$style(".butt{background:#0073b7;} 
                    .butt{color: white !important;}")))
  })
  
  #### Item to download Results #####
  output$downloadResults <- renderUI({
    if(is.null(rea.values$datasetProcess) || 
       length(rea.values$datasetProcess) != length(rea.values$datasetList)) return()
    
    column(
      width = 12, 
      downloadButton(outputId = "download", 
                     label = "Download results", class = "butt2"),
      tags$head(
        tags$style(".butt2{background:#0073b7;} 
                    .butt2{color: white !important;} ")))
  })
  
  
  
  #############################################
  # dynamic content #
  output$mycontent <- renderUI({
    
    items.list <- list()
    
    for(omics in c("RNAseq", "proteomics", "metabolomics")){
      
      items.list[[omics]] <-  lapply(1:10, function(i){
        tabItem(tabName = paste0(omics, "Analysis", i),
                uiOutput(paste0(omics, "AnalysisUI", i)))
      })
    }
    
    itemsOmics <- purrr::reduce(items.list, c)
    
    items <- c(
      
      list(
        
        #### Cover Page        ####
        ###########################
        tabItem(tabName = "coverPage",
                coverPageUI()
                
        ),
        #### Cover Page        ####
        ###########################
        tabItem(tabName = "GlossaryPage",
                GlossaryPageUI()
                
        ),
        #### Import data       ####
        ###########################
        tabItem(tabName = "importData",
                
                .modLoadDataUI("data")
        ),
        
        #### Set Up statistical model & hypothesis ####
        ###############################################
        tabItem(tabName = "SetUpModel",
                
                .modGLMmodelUI("model")
        )
      ),
      #### omics analysis ####
      ########################
      itemsOmics,
      list(
        
        #### analysis summary ####
        ########################
        tabItem(tabName = "omicsSum",
                uiOutput(outputId = "omicsSum_UI")
        ),
        
        #### MOFA ####
        ########################
        tabItem(tabName = "withMOFA",
                uiOutput(outputId = "withMOFA_UI")
        ),
        #### MixOmics ####
        ########################
        tabItem(tabName = "withMixOmics",
                # h5("in coming :)")
                uiOutput(outputId = "withMixOmics_UI") 
        )
      )
      
    )
    do.call(tabItems, items)
  })
  
  observe({
    
    lapply(names(rea.values$datasetList), function(omics){
      
      lapply(names(rea.values$datasetList[[omics]]), function(i){
        
        switch (omics,
                "RNAseq" = {output[[paste0("RNAseqAnalysisUI", i)]] <- renderUI({
                  
                  tabsetPanel(
                    
                    #### Data Exploratory & QC ####
                    ###############################
                    tabPanel("Pre-processing",
                             tags$br(),
                             tags$br(),
                             
                             QCNormalizationTabUI(paste0("RNAseq",i))
                             
                    ),
                    
                    #### Diff analysis  ####
                    ######################################
                    tabPanel("Differential analysis",
                             tags$br(),
                             tags$br(),
                             DiffExpAnalysisUI(paste0("RNAseq",i))
                    ),
                    #### Co-expression analysis  ####
                    ######################################
                    tabPanel("Co-expression analysis",
                             tags$br(),
                             tags$br(),
                             CoSeqAnalysisUI(paste0("RNAseq",i))
                             #verbatimTextOutput("Asuivre")
                    ),
                    #### enrichment analysis  CPR ####
                    ######################################
                    tabPanel("Annotation Enrichment",
                             tags$br(),
                             tags$br(),
                             .modEnrichmentUI(paste0("RNAseq",i))
                    )
                  )
                })},
                "proteomics" = {output[[paste0("proteomicsAnalysisUI", i)]] <- renderUI({
                  tabsetPanel(
                    
                    #### Data Exploratory & QC ####
                    ###############################
                    tabPanel("Pre-processing",
                             tags$br(),
                             tags$br(),
                             
                             QCNormalizationTabUI(paste0("proteomics",i))
                    ),
                    #### Diff analysis  ####
                    ######################################
                    tabPanel("Differential analysis",
                             tags$br(),
                             tags$br(),
                             DiffExpAnalysisUI(paste0("proteomics",i))
                    ),
                    #### Co-expression analysis  ####
                    ######################################
                    tabPanel("Co-expression analysis",
                             tags$br(),
                             tags$br(),
                             CoSeqAnalysisUI(paste0("proteomics",i))
                    ),
                    ### enrichment analysis CPR ####
                    #####################################
                    tabPanel("Annotation Enrichment",
                             tags$br(),
                             tags$br(),
                             .modEnrichmentUI(paste0("proteomics",i))
                    )
                  )
                })},
                "metabolomics" = {output[[paste0("metabolomicsAnalysisUI", i)]] <- renderUI({
                  
                  tabsetPanel(
                    
                    #### Data Exploratory & QC ####
                    ###############################
                    tabPanel("Pre-processing",
                             tags$br(),
                             tags$br(),
                             
                             QCNormalizationTabUI(paste0("metabolomics",i))
                    ),#### Diff analysis  ####
                    ######################################
                    tabPanel("Differential analysis",
                             tags$br(),
                             tags$br(),
                             DiffExpAnalysisUI(paste0("metabolomics",i))
                    ),
                    #### Co-expression analysis  ####
                    ######################################
                    tabPanel("Co-expression analysis",
                             tags$br(),
                             tags$br(),
                             CoSeqAnalysisUI(paste0("metabolomics",i))
                    ),
                    ### enrichment analysis CPR ####
                    #####################################
                    tabPanel("Annotation Enrichment",
                             tags$br(),
                             tags$br(),
                             .modEnrichmentUI(paste0("metabolomics",i))
                    )
                  )
                })},
        )
      })
    })
  })
  
  #### analysis Summary ####
  ###############################
  output$omicsSum_UI <- renderUI({
    
    .modSingleOmicAnalysesSummaryUI("omics")
  })

  #### MOFA data integration ####
  ###############################
  output$withMOFA_UI <- renderUI({
    
    .modIntegrationAnalysisUI("mofaSetting", method = "MOFA")
    
  })
  
  #### MixOmics data integration ####
  ###################################
  output$withMixOmics_UI <- renderUI({
    
    .modIntegrationAnalysisUI("mixomicsSetting", method = "mixOmics")
    
  })
  
  ########################################################################
  ######################### MAIN #########################################
  
  ##########################################
  # Part0 : presentation page
  ##########################################
  
  # # dir
  # shinyDirChoose(input = input, id = 'dir0', roots = c(home = '~'))
  # output$filepaths <- renderPrint({parseDirPath(roots = c(home = '~'), selection = input$dir0)})
  
  ##########################################
  # Part1 : load data
  ##########################################
  
  # load omics data and experimental design
  # set reference
  # set type of factor (bio/batch)
  # check design (complete and balanced)
  inputData <- callModule(.modLoadData, "data", rea.values)
  
  ##########################################
  # Part2 : Set GLM model
  ##########################################
  
  # display set up model Item
  # if no error message
  inputModel <- list()
  observeEvent(rea.values[["loadData"]], {
    
    #continue only if message is true or warning
    validate({
      need(rea.values$validate.status == 0, message = "set design step failed")
    })
    
    # display design menu
    output$SetUpModelMenu <- renderMenu({
      
      validate({
        need(rea.values$loadData == TRUE, message = FALSE)
      })
      menuItem(text = "Experimental Design", tabName = "SetUpModel", icon = icon('vials'))
    })
    
  }, ignoreInit = TRUE)
  
  # set GLM model
  # and select list of contrast to test
  {
    inputModel <- callModule(.modGLMmodel, "model", rea.values)
  }
  
  
  
  ##########################################
  # Part3 : ANALYSE
  ##########################################
  
  
  # for each omics data type
  # and for each dataser
  # if no error message
  observe({
    
    ##########################################
    # Part3 : Data Exploratory
    ##########################################
    lapply(names(rea.values$datasetList), function(omics){
      
      lapply(names(rea.values$datasetList[[omics]]), function(i){
        
        rea.values[[rea.values$datasetList[[omics]][[i]]]] <- reactiveValues(
          process    = FALSE,
          diffAnal   = FALSE,
          diffValid  = FALSE,
          coExpAnal  = FALSE,
          diffAnnot  = FALSE,
          coExpAnnot = FALSE,
          
          compCheck  = TRUE,
          message    = "",
          
          DiffValidContrast = NULL,
          CoExpClusterNames = NULL,
          omicsType = omics
        )
        
        ##########################################
        # Part3 : Data Exploratory
        ##########################################
        inputNorm <- callModule(module  = QCNormalizationTab, id = paste0(omics,i),
                                dataset = session$userData$FlomicsMultiAssay@metadata$omicList[[omics]][[i]],
                                rea.values = rea.values)
        
        ##########################################
        # Part5 :  Diff Analysis
        ##########################################
        inputDiff <- callModule(module  = DiffExpAnalysis, id = paste0(omics, i),
                                dataset = session$userData$FlomicsMultiAssay@metadata$omicList[[omics]][[i]],
                                rea.values = rea.values)
        
        ##########################################
        # Part6 : Co-Expression Analysis
        ##########################################
        callModule(module  = CoSeqAnalysis, id = paste0(omics, i),
                   dataset = session$userData$FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)
        
        ##########################################
        # Part7 : Enrichment Analysis CPR
        ##########################################
        callModule(module  = .modEnrichment, id = paste0(omics, i),
                   dataset = session$userData$FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)
        
      })
    })
    
  })
  
  callModule(module = .modSingleOmicAnalysesSummary, id = "omics", rea.values = rea.values)
  callModule(module = .modIntegrationAnalysis, id = "mixomicsSetting", rea.values = rea.values, method = "mixOmics")
  callModule(module = .modIntegrationAnalysis, id = "mofaSetting",     rea.values = rea.values, method = "MOFA")
 
  
  ##########################################
  # Part8 : RMD REPORT
  ##########################################
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"

    filename = function(){
      projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
      paste0(format(Sys.time(), "%Y_%m_%d"), "_", projectName, ".html")
    },
    content = function(file) {

      projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
      outDir <- file.path(tempdir(), paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName))
      dir.create(path = outDir, showWarnings=FALSE)
      
      generateReport(object = session$userData$FlomicsMultiAssay,
                     tmpDir = outDir, fileName = file)
    }
  )
  
  ##########################################
  # Part9 : Download results as an archive
  ##########################################

  output$download <- downloadHandler(
    # For PDF output, change this to "report.pdf"

    filename = function(){
      projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
      paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName, ".tar.gz")
    },
    content = function(file) {

      projectName  <- getProjectName(session$userData$FlomicsMultiAssay)
      outDir <- file.path(tempdir(), paste0(format(Sys.time(),"%Y_%m_%d"),"_", projectName))
      dir.create(path = outDir, showWarnings=FALSE)
      
      generateReport(object = session$userData$FlomicsMultiAssay,
                     tmpDir = outDir, archiveName = file, export = TRUE)
      
    })
  # # Automatically bookmark every time an input changes
  # observe({
  #   reactiveValuesToList(input)
  #   session$doBookmark()
  # })
  # # Update the query string
  # onBookmarked(updateQueryString)
}

