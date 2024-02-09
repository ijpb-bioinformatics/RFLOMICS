
##########################################
# Part3 : Data Exploratory & Pre-processing
##########################################

### UI
QCNormalizationTabUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = span(tagList(icon("filter"), "     Data Filtering and Normalization ", tags$small("(Scroll down for instructions)"))),
          width = 12, status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE)
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


### server
QCNormalizationTab <- function(input, output, session, dataset, rea.values){
  
  local.rea.values <- reactiveValues(message = NULL)
  
  #---- sample list----
  output$selectSamplesUI <- renderUI({
    
    sampleList <- colnames(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]])
    pickerInput(
      inputId  = session$ns("selectSamples"),
      label    = "Sample list:",
      choices  = sampleList,
      options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
      multiple = TRUE,
      selected = sampleList)
  })
  
  
  #---- Completeness----
  output$completenessUI <- renderUI({
    
    local.rea.values$message <- NULL
    completeCheckRes <- CheckExpDesignCompleteness(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], input$selectSamples)
    
    if(isTRUE(completeCheckRes[["error"]])){
      
      local.rea.values$message <- completeCheckRes[["messages"]]
    }
    # plot of count per condition
    if(!is.null(completeCheckRes[["plot"]])) {list(renderPlot(completeCheckRes[["plot"]]))}
    
  })
  
  #---- adapted parameters for each omics type ----
  output$paramUI <- renderUI({
    
    # setting for RNAseq
    paramRNAseq.list <- list(
      
      fluidRow(
        column(12, h5("Low count Filtering (CPM):")),
        column(8,
               selectInput(inputId  = session$ns("Filter_Strategy"),
                           label    = "Strategy",
                           choices  = c("NbConditions" = "NbConditions",  "NbReplicates" = "NbReplicates"),
                           selected = "NbReplicates")
        ),
        column(4,
               numericInput(inputId = session$ns("FilterSeuil"),
                            label="Cutoff:",
                            value=1, min = 1, max=10, step = 1 )
        )
      ),
      fluidRow(
        column(12,
               h5("Gene counts normalization:"),
               selectInput(inputId  = session$ns("selectNormMethod"),
                           label    = "method",
                           choices  =  list("TMM (edgeR)" = "TMM"),
                           selected = "TMM"))),
      fluidRow( column(12, actionButton(session$ns("run"),"Run")))
    )
    
    #---- setting for meta or RNAseq----
    paramProtMeta.list <- list(
      
      radioButtons(
        inputId  = session$ns("dataTransform"),
        label    = "Data transformation:",
        choices  = c("log1p" = "log1p","squareroot" = "squareroot",
                     "log2" = "log2", "log10" = "log10", "none" = "none"), 
        selected = "none"),
      hr(),
      
      radioButtons(inputId = session$ns("selectProtMetNormMethod"),
                   label    = "Normalization method:",
                   choices  =   c("median" = "median","totalSum" = "totalSum", "none" = "none"),
                   selected = "none"),
      
      
      actionButton(session$ns("run"),"Run")
    )
    
    switch(getOmicsTypes(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]]),
           "RNAseq"       = { paramRNAseq.list },
           "proteomics"   = { paramProtMeta.list },
           "metabolomics" = { paramProtMeta.list }
    )
  })
  
  #### PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "factors")
  
  #### PCA for metadata axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "meta")
  
  #---- dataset filtering summary----
  output$filtSummary1UI <- renderUI({
    
    SE.data  <- session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]]
    
    tagList(
      
      box(title = length(names(SE.data)), 
          width = 6, background = "maroon", 
          paste0("Initial number of ", omicsDic(SE.data)$variableName)
      ),
      box( title = length(colnames(SE.data)),
           width = 6, background = "light-blue", "Initial number of samples"
      )
    )
    
  })
  
  #---- Summary UI----
  output$filtSummary2UI <- renderUI({
    
    if (rea.values[[dataset]]$process == FALSE) return()
    
    SE.data  <- session$userData$FlomicsMultiAssay[[dataset]]
    
    tagList(
      box(
        title = length(names(SE.data)), width = 6, background = "fuchsia", 
        paste0("Number of filtered ", omicsDic(SE.data)$variableName)
      ),
      box(title = length(colnames(SE.data)),
          width = 6, background = "purple", "Number of filtered samples"
      )
    )
    
  })
  
  
  #---- Exploratory of Biological and Technical variability----
  output$tabPanelUI <- renderUI({
    
    if (rea.values$model == FALSE) return()
    
    MAE.data <- session$userData$FlomicsMultiAssay
    SE.data <- MAE.data[[paste0(dataset, ".raw")]]
    
    tabPanel.default.list <- list(
      
      tabPanel("Distribution (density)",
               tags$br(),
               uiOutput(session$ns("CountDistUI"))
      ),
      tabPanel("Distribution (boxplot)",
               tags$br(),
               uiOutput(session$ns("boxplotUI"))
      ),
      tabPanel("Principal component analysis",
               tags$br(),
               uiOutput(session$ns("PCAcoordUI"))
      )
    )
    
    tabPanel.list <- list()
    switch(getOmicsTypes(SE.data),
           "RNAseq" = {
             tabPanel.list <- c(
               list(
                 tabPanel("Library size",
                          tags$br(),
                          uiOutput(session$ns("LibSizeUI")))
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
  
  #---- QC plot for raw/processed data----
  # library size plot only for RNAseq data
  output$LibSizeUI <- renderUI({
    plot <- renderPlot(plotLibrarySize(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], raw = TRUE))
    
    if(rea.values[[dataset]]$process != FALSE){
      plot <- list(renderPlot(plotLibrarySize(session$userData$FlomicsMultiAssay[[dataset]])), plot)
    }
    
    return(plot)
  })
  
  # value (count/intensity) distribution (boxplot/density)
  output$boxplotUI <- renderUI({
    plot <- renderPlot(plotDataDistribution(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], plot = "boxplot", raw = TRUE))
    
    if(rea.values[[dataset]]$process != FALSE){
      plot <- list(renderPlot(plotDataDistribution(session$userData$FlomicsMultiAssay[[dataset]], plot = "boxplot", raw = FALSE)), plot)
    }
    
    return(plot)
  }) 
  
  output$CountDistUI <- renderUI({
    plot <- renderPlot(plotDataDistribution(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], plot = "density", raw = TRUE))
    
    if(rea.values[[dataset]]$process != FALSE){
      plot <- list(renderPlot(plotDataDistribution(session$userData$FlomicsMultiAssay[[dataset]], plot = "density", raw = FALSE)), plot)
    }
    
    return(plot)
  })   
  
  # PCA plot
  output$PCAcoordUI <- renderUI({
    
    factors_type   <- getFactorTypes(session$userData$FlomicsMultiAssay)
    choices        <- names(factors_type)
    names(choices) <- paste( names(factors_type),  paste0("(", factors_type, ")"))
    
    ui <- list(
      
      fluidRow(
        tags$br(),
        column(width = 6,
               radioButtons(inputId = session$ns("PCA.factor.condition"),
                            label = 'Levels:', choices = c("groups", choices), selected = "groups")),
        column(width = 6, UpdateRadioButtonsUI(session$ns("factors"))),
        tags$br(),
      ),
      plotOutput(session$ns("raw.PCAcoord")))
    
    if(rea.values[[dataset]]$process != FALSE){
      
      ui <- list(plotOutput(session$ns("norm.PCAcoord")), ui)
    }
    
    return(ui)
  })
  
  output$raw.PCAcoord <- renderPlot({
    
    PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
    PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
    condGroup <- input$PCA.factor.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], raw="raw", axes=c(PC1.value, PC2.value), groupColor=condGroup)
    
  })
  output$norm.PCAcoord <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
    PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
    condGroup <- input$PCA.factor.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[dataset]], raw = "norm", axes=c(PC1.value, PC2.value), groupColor=condGroup)
    
  })
  
  #---- run preprocessing - Normalization/transformation, filtering...----
  observeEvent(input$run, {
    
    # check if input$selectSamples is empty
    if(is.null(input$selectSamples)){
      showModal(modalDialog(title = "Error message", "please select samples."))
    }
    validate({ need(!is.null(input$selectSamples), message="please select samples.") })
    
    # check completeness for curent dataset
    if(!is.null(local.rea.values$message)){
      showModal(modalDialog(title = "Error message", local.rea.values$message))
    }
    validate({ need(is.null(local.rea.values$message), message=local.rea.values$message) })
    
    # get parameters 
    switch( getOmicsTypes(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]]),
            "RNAseq" = {
              param.list <- list(Filter_Strategy = input$Filter_Strategy, 
                                 CPM_Cutoff = input$FilterSeuil, 
                                 NormMethod=input$selectNormMethod)
            },
            {
              param.list <- list(transform_method = input$dataTransform,
                                 NormMethod=input$selectProtMetNormMethod)
            }
    )
    param.list <- c(param.list, list(samples = input$selectSamples))
    
    
    if(check_run_process_execution(session$userData$FlomicsMultiAssay, dataset = dataset, param.list = param.list) == FALSE &&
       rea.values[[dataset]]$process == TRUE) return()
    
    # re-initialize reactive values
    rea.values[[dataset]]$process   <- FALSE
    rea.values[[dataset]]$diffAnal  <- FALSE
    rea.values[[dataset]]$coExpAnal <- FALSE
    rea.values[[dataset]]$diffAnnot <- FALSE
    rea.values[[dataset]]$diffValid <- FALSE
    
    # re-initialize list of diff analized dataset
    if (dataset %in% rea.values$datasetDiff){
      rea.values$datasetDiff <- rea.values$datasetDiff[-which(rea.values$datasetDiff == dataset)]
    }
    if (dataset %in% rea.values$datasetProcess){
      rea.values$datasetProcess <- rea.values$datasetProcess[-which(rea.values$datasetProcess == dataset)]
    }
    
    # remove integration analysis
    session$userData$FlomicsMultiAssay <- resetFlomicsMultiAssay(session$userData$FlomicsMultiAssay, results = c("IntegrationAnalysis"))
    
    print(paste0("# 3  => Data processing: ", dataset))
   
    session$userData$FlomicsMultiAssay <- 
      runDataProcessing(object = session$userData$FlomicsMultiAssay, SE.name = dataset, samples = input$selectSamples,  
                        lowCountFiltering_strategy=param.list[["Filter_Strategy"]], lowCountFiltering_CPM_Cutoff=param.list[["CPM_Cutoff"]], 
                        normMethod=param.list[["NormMethod"]], transformMethod=param.list[["transform_method"]])
    
    rea.values[[dataset]]$process <- TRUE
    
    rea.values$datasetProcess <- unique(c(rea.values$datasetProcess , dataset))
    
  }, ignoreInit = TRUE)
  
  return(input)
  
}


############## functions ###############

# ----- check run norm execution ------
check_run_process_execution <- function(object.MAE, dataset, param.list = NULL){
  
  if(!dataset %in% names(object.MAE))  return(TRUE)
  
  SE <- object.MAE[[dataset]]
  
  if(!dplyr::setequal(SE$samples, param.list$samples)) return(TRUE)
  
  switch( getOmicsTypes(object.MAE[[dataset]]),
          "RNAseq" = {
            # filtering setting
            if(param.list$Filter_Strategy != getFilterSettings(SE)$filterStrategy) return(TRUE)
            if(param.list$CPM_Cutoff      != getFilterSettings(SE)$cpmCutoff)      return(TRUE)
            
            # normalisation setting
            if(param.list$NormMethod      != getNormSettings(SE)$method) return(TRUE)
          },
          "proteomics" = {
            
            if(param.list$transform_method != getTransSettings(SE)$method) return(TRUE)
            if(param.list$NormMethod       != getNormSettings(SE)$method) return(TRUE)
          },
          "metabolomics" = {
            
            if(param.list$transform_method != getTransSettings(SE)$method) return(TRUE)
            if(param.list$NormMethod       != getNormSettings(SE)$method) return(TRUE)
          }
  )
  
  return(FALSE)
}
