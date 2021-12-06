
library(shiny)
library(shinydashboard)
library(shinyFiles)

rm(list = ls())


shinyServer(function(input, output, session) {

    # This is to get the desired menuItem selected initially. 
    # selected=T seems not to work with a dynamic sidebarMenu.
    observeEvent(session, {
      updateTabItems(session = session, inputId = "tabs", selected = "coverPage")
    })
    
    # reactive value for reinitialisation of UIoutput
    rea.values <- reactiveValues(
      loadData = FALSE,
      model    = FALSE,
      analysis = FALSE
    )
    
    # Use reactive values for dynamic Items
    #subitems <- reactiveVal(value = NULL)
    
    # dynamic sidebar menu #
    output$mysidebar <- renderUI({
    
      tagList(
          sidebarMenu(id="tabs",
                      menuItem(text = "Welcome", tabName = "coverPage", icon = icon('dna'), selected = TRUE),
                      menuItem(text = "Load Data", tabName = "importData", icon = icon('download')),
                      menuItemOutput(outputId = "SetUpModelMenu"),
                      menuItemOutput(outputId = "omics"),
                      menuItemOutput(outputId = "Integration")
          ),
          
          tags$br(),
          tags$br(),
          uiOutput("runReport")
        
          #downloadButton(outputId = "report", label = "Generate report")
          )
    })
    
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
        
        itemsOmics,
        list(
          tabItem(tabName = "OmicsIntegration",
                  h5("in coming :)")
          )
        )
        
      )
      do.call(tabItems, items)
    })
    
    lapply(1:10, function(i){
      
      output[[paste0("RNAseqAnalysisUI", i)]] <- renderUI({
        
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
      })
    })
    lapply(1:10, function(j){  
      
      output[[paste0("proteomicsAnalysisUI", j)]] <- renderUI({
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
      })
    })
    lapply(1:10, function(k){
      
      output[[paste0("metabolomicsAnalysisUI", k)]] <- renderUI({
        
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
      })
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
    inputData <- callModule(LoadOmicsData, "data", rea.values)
    

    ##########################################
    # Part2 : Set GLM model
    ##########################################

    # display set up model Item
    # if no error message
    inputModel <- list()
    observeEvent(inputData[["loadData"]], {
      
      #subitems(NULL)
      
      #continue only if message is true or warning
      validate({
        need(rea.values$validate.status == 0, message="set design step failed")
      })

      # display design menu
      output$SetUpModelMenu <- renderMenu({
        
        menuItem(text = "Experimental Design", tabName = "SetUpModel", icon = icon('vials'))
      })

      
    })

    # set GLM model
    # and select list of contrast to test
    inputModel <- callModule(GLM_model, "model", rea.values)

    ##########################################
    # Part3 : ANALYSE
    ##########################################
    
    #### Item for each omics #####
    # display omics Item
    output$omics <- renderMenu({
      
      if(rea.values$analysis == FALSE) return()
      
      menu_list <- list()
      menu_list <- list(
        menu_list,
        sidebarMenu(id = "sbm",
                    
                    #lapply(subitems(), function(omics){
                    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
                      
                      do.call(what = menuItem,
                              args = c(text = paste0(omics, " Analysis"), tabName = paste0(omics, "Analysis"), icon = icon('chart-line'),
                                       
                                       lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                                         menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]),
                                                     tabName = paste0(omics, "Analysis", i), icon = icon('chart-area'))
                                       })
                              )
                      )
                    })
        )
      )
      sidebarMenu(.list = menu_list)
    })
    # update Item menu
    #subitems(names(FlomicsMultiAssay@metadata$omicList))
    #updateTabItems(session, "sbm", selected = "SetUpModelMenu")
    
    
    output$runReport <- renderUI({
      if(rea.values$analysis == FALSE) return()
      
      downloadButton(outputId = "report", label = "Generate report")
    })
    
    
    output$Integration <- renderMenu({
      if(rea.values$analysis == FALSE) return()
      
      menuItem(text = "Data Integration", tabName = "OmicsIntegration", icon = icon('network-wired'))  
    })
    
    # for each omics data type
    # and for each dataser
    # if no error message
    
    observeEvent(inputModel$validContrasts, {
      
      
      ##########################################
      # Part3 : Data Exploratory
      ##########################################
      #inputNorm <- list()
      lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
        
        lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
          
          rea.values[[FlomicsMultiAssay@metadata$omicList[[omics]][[i]]]] <- reactiveValues(
            diffAnal   = FALSE,
            coExpAnal  = FALSE,
            diffAnnot  = FALSE,
            diffValid  = FALSE,
            coExpAnnot = FALSE
          )
          inputNorm <- callModule(module  = QCNormalizationTab, id = paste0(omics,i),
                                                      dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], 
                                                      rea.values = rea.values)
          inputDiff <- callModule(module  = DiffExpAnalysis, id = paste0(omics, i), 
                                                       dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], 
                                                       rea.values = rea.values)
          
          callModule(module  = CoSeqAnalysis, id = paste0(omics, i),
                     dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)
      
      
          callModule(module  = AnnotationEnrichment, id = paste0(omics, i),
                     dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)
          
          
        })
      })
      
    }, ignoreInit = TRUE)
    
    
    
    
  #   observeEvent(inputModel$validContrasts, {
  # 
  # 
  #   ##########################################
  #   # Part3 : Data Exploratory
  #   ##########################################
  #   inputNorm <- list()
  #   lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
  # 
  #     lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
  #       
  #       #rea.values[[FlomicsMultiAssay@metadata$omicList[[omics]][[i]]]] <- reactiveVal()
  #       inputNorm[[paste0(omics, i)]] <- callModule(module  = QCNormalizationTab, id = paste0(omics,i),
  #                                                   dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], 
  #                                                   rea.values = rea.values)
  #       })
  #   })
  # 
  # #   ##########################################
  # #   # Part5 :  Diff Analysis
  # #   ##########################################
  # # 
  # #   inputDiff <- list()
  # #   lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
  # # 
  # #     lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
  # # 
  # #       #observeEvent(inputModel$validContrasts, {
  # # 
  # #           inputDiff[[paste0(omics, i)]] <<- callModule(module  = DiffExpAnalysis, id = paste0(omics, i), 
  # #                                                        dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], 
  # #                                                        rea.values = rea.values)
  # # 
  # #       #})
  # #     })
  # #   })
  # # 
  # # 
  # # ##########################################
  # # # Part6 : Co-Expression Analysis
  # # ##########################################
  # #   lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
  # # 
  # #     lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
  # # 
  # #       observeEvent(inputDiff[[paste0(omics, i)]]$validContrast, {
  # # 
  # #         callModule(module = CoSeqAnalysis, id = paste0(omics, i),
  # #                    dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)
  # #       }, ignoreInit = TRUE)
  # #     })
  # #   })
  # # 
  # # 
  # #   ##########################################
  # #   # Part7 : Enrichment Analysis
  # #   ##########################################
  # #   lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
  # # 
  # #     lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
  # # 
  # #       observeEvent(inputDiff[[paste0(omics, i)]]$validContrast, {
  # # 
  # #         callModule(module = AnnotationEnrichment, id = paste0(omics, i),
  # #                    dataset = FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)
  # #       })
  # #     })
  #   })
 
    
    

    ##########################################
    # Part8 : RMD REPORT
    ##########################################

    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = paste0(FlomicsMultiAssay@metadata$projectName, "_", format(Sys.time(), "%Y_%m_%d_%H_%M"), ".html"),
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).

        tempReport <-  "report.Rmd" # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #tempReport <- file.path(tempdir(), "report.Rmd")
        #file.copy("report.Rmd", tempReport, overwrite = TRUE)

        # TEST
        # save FE object in .Rdata and load it during report execution
        save(FlomicsMultiAssay,file=file.path(tempdir(), "FlomicsMultiAssay.RData"))

        # Set up parameters to pass to Rmd document
        params <- list( FEdata = file.path(tempdir(), "FlomicsMultiAssay.RData"),
                        pngDir = tempdir())

        print(tempdir())
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv()))

      }
    )

    # # Automatically bookmark every time an input changes
    # observe({
    #   reactiveValuesToList(input)
    #   session$doBookmark()
    # })
    # # Update the query string
    # onBookmarked(updateQueryString)
 })
 
