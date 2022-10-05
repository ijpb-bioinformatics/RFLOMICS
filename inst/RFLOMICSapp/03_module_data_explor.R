
##########################################
# Part3 : Data Exploratory & Pre-processing
##########################################


QCNormalizationTabUI <- function(id){

  #name space for id
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = span(tagList(icon("filter"), "     Data Filtering and Normalization ", tags$small("(Scroll down for instructions)"))),
          width = 12, status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,

          p("For each diagnostic plot, both raw and processed (filtered, normalized, ...) data are displayed with expertised default parameters."),
          p("- You may first have a look at to default processed plots to eventually identify outliers sample"),
          p("- You may check onto the ACP factorial map to control that samples group with the replicat of the same condition."),
          p("- You may play with filtered parameters to see if it improved the grouping"),
          p("- If not, you may remove outliers from the sample list, update the analysis and check again."),
          p("- To quickly overview the % of variability which is associated to each design factors, go to PCA (2/2) ")
      )
    ),
    fluidRow(
      column(4,
             box(width = 14, title = span(tagList(icon("sliders-h"), "  ", "Setting")), status = "warning",
                 uiOutput(ns("selectSamplesUI")),
                 uiOutput(ns("completenessUI")),
                 hr(),
                 uiOutput(ns("paramUI"))),
             fluidRow(
               uiOutput(ns("filtSummary1UI"))),
             fluidRow(
               uiOutput(ns("filtSummary2UI")))),
      
      column(8, uiOutput(ns("tabPanelUI")))
    )
  )
}



QCNormalizationTab <- function(input, output, session, dataset, rea.values){

  local.rea.values <- reactiveValues(dataset.processed.SE = NULL, 
                                     dataset.raw.SE = NULL,
                                     compCheck = TRUE,
                                     message = "")
  
  ### check completeness
  completeCheckRes <- CheckExpDesignCompleteness(object = session$userData$FlomicsMultiAssay,
                                                 sampleList = colnames(session$userData$FlomicsMultiAssay[[dataset]]))
  # stock message in MAE
  session$userData$FlomicsMultiAssay@metadata$completeCheck[["error"]]   <- completeCheckRes[["error"]]
  session$userData$FlomicsMultiAssay@metadata$completeCheck[["warning"]] <- completeCheckRes[["warning"]]
  # reactive values
  rea.values[[dataset]]$compCheck <- TRUE
  if(!is.null(completeCheckRes[["error"]])){ rea.values[[dataset]]$compCheck <- FALSE }
  
  
  #### sample list  ####
  output$selectSamplesUI <- renderUI({

    sampleList <- colnames(session$userData$FlomicsMultiAssay[[dataset]])
      # sample list :
      pickerInput(
        inputId  = session$ns("selectSamples"),
        label    = "Sample list :",
        choices  = sampleList,
        options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
        multiple = TRUE,
        selected = sampleList)
  })
  
  #### check completeness ####
  output$completenessUI <- renderUI({
    
    completeCheckRes <- CheckExpDesignCompleteness(session$userData$FlomicsMultiAssay, input$selectSamples)
    
    local.rea.values$compCheck <- TRUE
    
    if(!is.null(completeCheckRes[["error"]])){
      local.rea.values$compCheck <- FALSE
      local.rea.values$message   <- completeCheckRes[["error"]]
      
      showModal(modalDialog(title = "Error message", completeCheckRes[["error"]]))
    }
    list(
         # plot of count per condition
         renderPlot( completeCheckRes[["plot"]])
    )
  })

  #### adapted parameters for each omics type
  output$paramUI <- renderUI({

    # setting for RNAseq
    paramRNAseq.list <- list(
      
      fluidRow(
        column(12, h5("Low count Filtering (CPM) :")),
        column(8,
               selectInput(inputId  = session$ns("Filter_Strategy"),
                           label    = "Stategy",
                           choices  = c("NbConditions" = "NbConditions",  "NbReplicates" = "NbReplicates"),
                           selected = "NbConditions")
               ),
        column(4,
               numericInput(inputId = session$ns("FilterSeuil"),
                            label="Cutoff :",
                            value=5, min = 1, max=10, step = 1 )
               )
      ),
      fluidRow(
          column(12,
             h5("Gene count normalization :"),
             selectInput(inputId  = session$ns("selectNormMethod"),
                         label    = "method",
                         choices  =  list("TMM (edgeR)" = "TMM"),
                         selected = "TMM"))),
      fluidRow( column(12, actionButton(session$ns("Update"),"Update")))
    )

    # setting for meta or RNAseq
    paramProtMeta.list <- list(
        # radioButtons(
        #   inputId  =session$ns("dataImputation"),
        #   label = "does the data need to be impute (no implemented yet) ?",
        #   choices=c("yes"="yes","no"="no"),
        #   selected="no"),
        # hr(),

        radioButtons(
          inputId  =session$ns("dataTransform"),
          label = "Does the data need to be transformed ?",
          choices=c("log1p"="log1p","squareroot"="squareroot","no"="none"),
          selected="none"),
        hr(),

        actionButton(session$ns("Update"),"Update")
    )

    switch (session$userData$FlomicsMultiAssay[[dataset]]@metadata$omicType,
            "RNAseq"       = { paramRNAseq.list },
            "proteomics"   = { paramProtMeta.list },
            "metabolomics" = { paramProtMeta.list }
    )
  })
  
  #### PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "rawData")
  # select factors for color PCA plot
  callModule(RadioButtonsCondition, "rawData")

  #### dataset filterin summary
  output$filtSummary1UI <- renderUI({
    
    tagList(
      
      box(title = length(names(session$userData$FlomicsMultiAssay[[dataset]])), 
          width = 6, background = "maroon", "nbr of entities before filtering"
      ),
      box( title = length(names(local.rea.values$dataset.processed.SE)), 
           width = 6, background = "light-blue", "nbr of entities after filtering"
      )
    )
    
  })
  output$filtSummary2UI <- renderUI({
    
    if(rea.values[[dataset]]$omicsType != "RNAseq") return()
    tagList(
      
      box(
        title = length(session$userData$FlomicsMultiAssay[[dataset]]@metadata$rowSums.zero), 
        width = 6, background = "fuchsia", "nbr of unexpressed genes"
      ),
      box(title = length(local.rea.values$dataset.processed.SE@metadata$FilteredFeatures), 
          width = 6, background = "purple", "nbr of low expressed genes"
      )
    )
    
  })
  #### processing data with default param
  # RNAseq : normalization
  # prot/meta : transformation
  param.list <- switch( session$userData$FlomicsMultiAssay[[dataset]]@metadata$omicType,
              "RNAseq"       = { list(Filter_Strategy = "NbConditions", CPM_Cutoff = 5, NormMethod = "TMM") },
              "proteomics"   = { list(transform_method = "none") },
              "metabolomics" = { list(transform_method = "none") } )

  print("# 5- Compute PCA for raw counts")
  session$userData$FlomicsMultiAssay[[dataset]] <-  RFLOMICS::RunPCA(session$userData$FlomicsMultiAssay[[dataset]])
  session$userData$FlomicsMultiAssay[[dataset]] <- session$userData$FlomicsMultiAssay[[dataset]]
  
  # run pre-preprocessing
  print(paste0("# 6- Data processing... ", dataset))
  SE.processed <- process_data(SE = session$userData$FlomicsMultiAssay[[dataset]], dataset = dataset,
                               samples = colnames(session$userData$FlomicsMultiAssay[[dataset]]), param.list = param.list)
  session$userData$FlomicsMultiAssay@ExperimentList[[paste0(dataset, ".filtred")]] <- SE.processed
  local.rea.values$dataset.processed.SE <- SE.processed
  
  SE.name <- paste0(dataset,".filtred")

  # remove SE processed if existe
  if (SE.name %in% names(session$userData$FlomicsMultiAssay)){
    session$userData$FlomicsMultiAssay <- session$userData$FlomicsMultiAssay[,, -which(names(session$userData$FlomicsMultiAssay) == SE.name)]
  }

  session$userData$FlomicsMultiAssay <- eval(parse(text = paste0('c( session$userData$FlomicsMultiAssay ,', SE.name, ' = SE.processed )')))

  
  #### Exploratory of Biological and Technical variability
  output$tabPanelUI <- renderUI({

    if(rea.values$model == FALSE) return()


    #if(is.null(FlomicsMultiAssay)) return()

    tabPanel.default.list <- list(

        tabPanel("Distribution (density)",
                  tags$br(),
                  tags$br(),
                  column(width = 12, plotOutput(session$ns("raw.CountDist"))),
                  column(width = 12, plotOutput(session$ns("norm.CountDist")))
        ),
        tabPanel("Distribution (boxplot)",
                 tags$br(),
                 tags$br(),
                 column(width = 12, plotOutput(session$ns("raw.boxplot"))),
                 column(width = 12, plotOutput(session$ns("norm.boxplot")))
        ),
        tabPanel("Principal component analysis",
                 tags$br(),
                 column(width = 12,  plotOutput(session$ns("raw.PCAcoord"))),
                 hr(),
                 fluidRow(
                   column(width = 6, RadioButtonsConditionUI(session$ns("rawData"))),
                   column(width = 6, UpdateRadioButtonsUI(session$ns("rawData")))
                 ),
                 hr(),
                 column(width = 12, plotOutput(session$ns("norm.PCAcoord")))#,
                 # tags$br(),
                 # column(width = 12, actionButton(session$ns("screenshotPCA_QC"),"Screenshot"))
        )
        # ,
        # tabPanel("Principal component analysis (2/2)",
        #          plotOutput(session$ns("QCdesignPCA"))
        #          )
    )

    tabPanel.list <- list()
     switch(session$userData$FlomicsMultiAssay[[dataset]]@metadata$omicType,
            "RNAseq" = {
              tabPanel.list <- c(
                list(
                  tabPanel("Library size",
                         # library size plot
                         #column(4, plotOutput(ns("LibSize"),   height = "400%")),
                         column(12, plotOutput(session$ns("raw.LibSize"))),
                         column(12, plotOutput(session$ns("norm.LibSize"))))
                  ),
                  tabPanel.default.list)
            },
            "proteomics" = {
              tabPanel.list <- tabPanel.default.list
            },
            "metabolomics" = {
              tabPanel.list <- tabPanel.default.list
            }
    )

    # Exploratory of Biological and Technical variability
    box(width = 14, title = "Exploratory of Biological and Technical variability", 
        solidHeader = TRUE, status = "warning",
      do.call(what = tabsetPanel, args = tabPanel.list))
  })

  #### QC plot for raw data
  # library size plot only for RNAseq data
  output$raw.LibSize <- renderPlot({

    Library_size_barplot.plot(session$userData$FlomicsMultiAssay[[dataset]])
  })
  
  # value (count/intensity) distribution (boxplot/density)
  output$raw.boxplot <- renderPlot({
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[dataset]], plot = "boxplot")
  })
  output$raw.CountDist <- renderPlot({
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[dataset]], plot = "density")
  })
  
  # PCA plot
  output$raw.PCAcoord <- renderPlot({
    #session$userData$FlomicsMultiAssay[[dataset]] <-  RunPCA(session$userData$FlomicsMultiAssay[[dataset]])
    PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    condGroup <- input$`rawData-condColorSelect`[1]
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[dataset]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
    
  })
  
  #### update pre-processing with new param
  observeEvent(input$Update, {

    #print(paste0("Update ",input$Update))

    # continue only if message is true or warning
    #completeCheckRes <- CheckExpDesignCompleteness(session$userData$FlomicsMultiAssay, input$selectSamples)
    
    if(isFALSE(local.rea.values$compCheck)){
      showModal(modalDialog(title = "Error message", local.rea.values$message))
    }
    validate({ need(!isFALSE(local.rea.values$compCheck), message=local.rea.values$message) })
    
    rea.values[[dataset]]$compCheck <- TRUE
    
    # re-initialize reactive values
    rea.values[[dataset]]$diffAnal  <- FALSE
    rea.values[[dataset]]$coExpAnal <- FALSE
    rea.values[[dataset]]$diffAnnot <- FALSE
    rea.values[[dataset]]$diffValid <- FALSE

    print("# 6.bis => Update data processing...")
    # re-initialize MAE object
    if (dataset %in% rea.values$datasetDiff){
      rea.values$datasetDiff <- rea.values$datasetDiff[-which(rea.values$datasetDiff == dataset)]
      }

    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]]@metadata$DiffExpAnal       <- list()
    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]]@metadata$CoExpAnal         <- list()
    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]]@metadata$DiffExpEnrichAnal <- list()
    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]]@metadata$CoExpEnrichAnal   <- list()

    # new selected param
    param.list <- list()
    switch( rea.values[[dataset]]$omicsType,
            "RNAseq" = {
              param.list <- list(Filter_Strategy = input$Filter_Strategy, CPM_Cutoff = input$FilterSeuil, NormMethod=input$selectNormMethod)
            },
            "proteomics" = {
              param.list <- list(transform_method = input$dataTransform)
            },
            "metabolomics" = {
              param.list <- list(transform_method = input$dataTransform)
            }
    )
    
    print(paste0("# 6- Data processing... ", dataset))
    processed.SE <- process_data(SE = session$userData$FlomicsMultiAssay[[dataset]], dataset = dataset, samples = input$selectSamples, param.list = param.list)

    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]] <- processed.SE
    
    local.rea.values$dataset.processed.SE <- processed.SE
    
    ## add new SE with processed data
    SE.name <- paste0(dataset,".filtred")
    
    # remove SE processed if existe
    if (SE.name %in% names(session$userData$FlomicsMultiAssay)){
      session$userData$FlomicsMultiAssay <- session$userData$FlomicsMultiAssay[,, -which(names(session$userData$FlomicsMultiAssay) == SE.name)]
    }
    
    session$userData$FlomicsMultiAssay <- eval(parse(text = paste0('c( session$userData$FlomicsMultiAssay ,', SE.name, ' = processed.SE )')))
    
  }, ignoreInit = TRUE)

  
  #### QC plot for processed data
  # library size plot only for RNAseq data
  output$norm.LibSize <- renderPlot({
    
    if(rea.values[[dataset]]$omicsType != "RNAseq") return()
    
    Library_size_barplot.plot(local.rea.values$dataset.processed.SE)
  })
  
  # value (count/intensity) distribution (boxplot/density)
  output$norm.boxplot <- renderPlot({
    if(rea.values[[dataset]]$omicsType != "RNAseq" && local.rea.values$dataset.processed.SE@metadata$transform_method == "none") return()
    
    Data_Distribution_plot(local.rea.values$dataset.processed.SE, plot = "boxplot")
  })
  
  output$norm.CountDist <- renderPlot({
    if(rea.values[[dataset]]$omicsType != "RNAseq" && local.rea.values$dataset.processed.SE@metadata$transform_method == "none") return()
    
    Data_Distribution_plot(local.rea.values$dataset.processed.SE, plot = "density")
  })
  
  # PCA plot
  output$norm.PCAcoord <- renderPlot({
    
    if(rea.values[[dataset]]$omicsType != "RNAseq" && local.rea.values$dataset.processed.SE@metadata$transform_method == "none") return()
    #session$userData$FlomicsMultiAssay[[dataset]] <-  RunPCA(session$userData$FlomicsMultiAssay[[dataset]])
    PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    condGroup <- input$`rawData-condColorSelect`[1]

    RFLOMICS::plotPCA(local.rea.values$dataset.processed.SE, PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup)

  })
  

    # 
    # # => display only for normalized RNAseq data
    # if (local.rea.values$dataset.processed.SE@metadata$omicType == "RNAseq"){
    #   output$norm.PCAcoord <- renderPlot({
    #     #local.rea.values$dataset.processed.SE <-  RunPCA(local.rea.values$dataset.processed.SE)
    #     PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    #     PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    #     condGroup <- input$`rawData-condColorSelect`[1]
    #     
    #     RFLOMICS::plotPCA(local.rea.values$dataset.processed.SE, PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
    #   })
    # }
    # 
    # 
    # # # save current PCA plot with fixed axis & color
    # # ## screenShot
    # # observeEvent(input$screenshotPCA_QC, {
    # #
    # #   PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    # #   PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    # #   condGroup <- input$`rawData-condColorSelect`[1]
    # #
    # #   # file.copy(file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
    # #   #           file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
    # #   #           overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
    # # })
    # 
    # #### PCA analysis QCdesign ####
    # output$QCdesignPCA <- renderPlot({
    #   
    #   # if(is.null(session$userData$FlomicsMultiAssay[[dataset]][["PCAlist"]][["raw"]])){
    #   #   session$userData$FlomicsMultiAssay[[dataset]] <-  RunPCA(session$userData$FlomicsMultiAssay[[dataset]])
    #   # }
    #   
    #   mvQCdesign(session$userData$FlomicsMultiAssay,data=dataset,PCA="raw", axis=5,
    #              pngFile=file.path(tempdir(), paste0(dataset,"_PCAdesignCoordRaw.png")))
    # })
    
    # #### PCA analysis QCdata ####
    # output$QCdata <- renderPlot({
    #
    #   mvQCdata(session$userData$FlomicsMultiAssay,data=dataset,PCA="raw",axis=5,
    #            pngFile=file.path(tempdir(), paste0(dataset,"_PCAmetaCorrRaw.png")))
    # })
    
  # # save current PCA plot with fixed axis & color
  # ## screenShot
  # observeEvent(input$screenshotPCA_Norm, {
  #   PC1.value <- as.numeric(input$`normData-Firstaxis`[1])
  #   PC2.value <- as.numeric(input$`normData-Secondaxis`[1])
  #   condGroup <- input$`normData-condColorSelect`[1]
  #
  #   file.copy(file.path(tempdir(), paste0(dataset,".filtred","_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
  #             file.path(tempdir(), paste0(dataset,".filtred","_PCAdesign_norm_PC",     PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
  #             overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
  # })
  #
  #
  #


    return(input)

}


############## functions ###############
###
process_data <- function(SE, dataset, samples , param.list = list(Filter_Strategy = "NbConditions", CPM_Cutoff = 1, NormMethod = "TMM", transform_method = "none")){

  print("# => select samples")
  SE.new <- SE[, SE$primary %in% samples]

  SE.new@metadata$Groups <- dplyr::filter(SE@metadata$Groups, samples %in% SE.new$primary)

   switch(SE.new@metadata$omicType,
          "RNAseq" = {
             #### Filter low abundance ####
             print("# => Low counts Filtering...")
             SE.processed <- FilterLowAbundance(SE.new, param.list[["Filter_Strategy"]], param.list[["CPM_Cutoff"]])

             #### Run Normalisation ####
             print("# => Counts normalization...")
             SE.processed <- RunNormalization(SE.processed, param.list[["NormMethod"]])
          },
         "proteomics" = {
             print("# => transformation data...")
             SE.processed <- TransformData(SE.new, transform_method = param.list[["transform_method"]])
             
         },
         "metabolomics" = {
             print("# => transformation data...")
             SE.processed <- TransformData(SE.new, transform_method = param.list[["transform_method"]])
         }
       )

  #### Run PCA for filtred & normalized data ####
  print("# => Compute PCA ")
  SE.processed <- RunPCA(SE.processed)
  
  return(SE.processed)
}


