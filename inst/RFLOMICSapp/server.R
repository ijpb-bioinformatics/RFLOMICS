
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
      analysis = FALSE,

      datasetList  = NULL,
      contrastList = NULL,
      datasetDiff  = NULL
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
                  h5("in coming :)")
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
                #### enrichment analysis  ####
                ######################################
                tabPanel("Annotation Enrichment",
                         tags$br(),
                         tags$br(),
                         AnnotationEnrichmentUI(paste0("RNAseq",i))
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
                         #verbatimTextOutput("Asuivre")
                ),
                #### enrichment analysis  ####
                ######################################
                tabPanel("Annotation Enrichment",
                         tags$br(),
                         tags$br(),
                         AnnotationEnrichmentUI(paste0("proteomics",i))
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
                         #verbatimTextOutput("Asuivre")
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

      omics_data_analysis_summaryUI("omics")

    })


    #### MOFA data integration ####
    ###############################
    output$withMOFA_UI <- renderUI({

      MOFA_settingUI("mofaSetting")

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

        validate({
          need(rea.values$loadData == TRUE, message="")
        })
        menuItem(text = "Experimental Design", tabName = "SetUpModel", icon = icon('vials'))
      })

    }, ignoreInit = TRUE)

    # set GLM model
    # and select list of contrast to test
    observe({
      inputModel <- callModule(GLM_model, "model", rea.values)
    })



    ##########################################
    # Part3 : ANALYSE
    ##########################################

    #### Item for each omics #####
    # display omics Item
    output$omics <- renderMenu({

      validate({
        need(rea.values$analysis == TRUE, message="")
      })

      menuItem(text = "Omics Analysis", tabName = "OmicsAnalysis", icon = icon('chart-area'),
               #icon = icon('chart-line'),
           lapply(names(rea.values$datasetList), function(omics){

             lapply(names(rea.values$datasetList[[omics]]), function(i){
                menuSubItem(text = rea.values$datasetList[[omics]][[i]],
                            tabName = paste0(omics, "Analysis", i))
              })
           })
      )




    })

    # observe({
    #
    #   updateTabItems(session, "tabs", selected = paste0(names(FlomicsMultiAssay@metadata$omicList)[1],
    #                                                     "Analysis",
    #                                                     names(FlomicsMultiAssay@metadata$omicList[[1]])[1]))
    # })


    # update Item menu
    #subitems(names(FlomicsMultiAssay@metadata$omicList))
    #updateTabItems(session, "sbm", selected = "SetUpModelMenu")


    #### Item for each data integration tools #####
    # display tool Item
    output$Integration <- renderMenu({

      validate({
        need(rea.values$analysis == TRUE && length(rea.values$datasetDiff) >=2, message = "")
      })


      menuItem(text = "Data Integration", tabName = "OmicsIntegration", icon = icon('network-wired'), startExpanded = FALSE,selected = FALSE,
           menuSubItem(text = "Dataset analysis summary", tabName = "omicsSum" ),
           menuSubItem(text = "with MOFA", tabName = "withMOFA" ),
           menuSubItem(text = "with MixOmics", tabName = "withMixOmics")
      )
    })

    #### Item for report #####
    output$runReport <- renderUI({
      if(rea.values$analysis == FALSE) return()

      downloadButton(outputId = "report", label = "Generate report")
    })



    # for each omics data type
    # and for each dataser
    # if no error message
    observe({

      if(rea.values$analysis == FALSE) return()
      ##########################################
      # Part3 : Data Exploratory
      ##########################################
      lapply(names(rea.values$datasetList), function(omics){

        lapply(names(rea.values$datasetList[[omics]]), function(i){

          rea.values[[rea.values$datasetList[[omics]][[i]]]] <- reactiveValues(
            diffAnal   = FALSE,
            diffValid  = FALSE,
            coExpAnal  = FALSE,
            diffAnnot  = FALSE,
            coExpAnnot = FALSE,

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
          # Part7 : Enrichment Analysis
          ##########################################
          callModule(module  = AnnotationEnrichment, id = paste0(omics, i),
                     dataset = session$userData$FlomicsMultiAssay@metadata$omicList[[omics]][[i]], rea.values = rea.values)


        })
      })


      callModule(module = omics_data_analysis_summary, id = "omics", rea.values = rea.values)

      callModule(module = MOFA_setting, id = "mofaSetting", rea.values = rea.values)

    })



    ##########################################
    # Part8 : RMD REPORT
    ##########################################

    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      
      filename = function(){
        projectName <- session$userData$FlomicsMultiAssay@metadata$projectName
        paste0(projectName, "_", format(Sys.time(), "%Y_%m_%d_%H_%M"), ".html")
        },
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        print(paste0("# ??- Creat html report... ", file))
        tempReport <-  paste0(path.package("RFLOMICS"), "/RFLOMICSapp/","report.Rmd") # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #tempReport <- file.path(tempdir(), "report.Rmd")
        #file.copy("report.Rmd", tempReport, overwrite = TRUE)

        # TEST
        # save FE object in .Rdata and load it during report execution
        projectName  <- session$userData$FlomicsMultiAssay@metadata$projectName
        rflomics.MAE <- session$userData$FlomicsMultiAssay
        RData.name   <- paste0(projectName, ".MAE.RData")
        save(rflomics.MAE, file=file.path(tempdir(), RData.name))

        # Set up parameters to pass to Rmd document
        print(file.path(tempdir(), RData.name))
        params <- list( FEdata = file.path(tempdir(), RData.name),
                        title  = paste0(projectName, "project"),
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

