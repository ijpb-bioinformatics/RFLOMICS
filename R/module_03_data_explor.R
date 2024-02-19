
##########################################
# Part3 : Data Exploratory & Pre-processing
##########################################

### UI
QCNormalizationTabUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = span(tagList(icon("filter"), "     Data exploration and pre-processing ", tags$small("(Scroll down for instructions)"))),
          width = 12, status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          p("Very important and non trivial step which are crucial for single-omics analysis and also for the integration of different omics data. Default settings have been expertized"),
          h4(tags$span("Explore:", style = "color:orange")),
          p("Identify the noise and the source of technical and biological variability thanks to the",tags$b("PCA"),""),
          h4(tags$span("Filter:", style = "color:orange")),
          p("Remove outliers samples so as low expressed entities which introduce noise in the data"),
          h4(tags$span("Transform and normalize:", style = "color:orange")),
          p("Transform data to linearize and make it more Gaussian-like and Normalize to identify and correct technical biases and make the data comparable across samples. Depending on the omics type, the pre-processing steps will be different: "),
          h5(tags$span("RNAseq data:", style = "color:blue")),
          p("",tags$b("Details on filtering step:"),""),
          p("By default, non expressed/non detected genes are removed"),
          p("By default, low expressed genes are removed according to their CPM. By default, genes with a cpm >= 1 in at least n Samples = min(NbReplicates) samples are kept."),
          p("NB: you can choose the other strategy which is to remove genes according to a cpm >= 1 in at least n Samples = NbConditions cpm threshold could be changed"),
          p("",tags$b("Details on normalization:"),""),
          p("TMM from edgeR is the default method. It was found to be the best methods (Dillies et al., 2013)"),
          h5(tags$span("Proteomic and metabolomic data", style = "color:blue")),
          p("",tags$b("Details on transformation:"),""),
          p("Log2 is the default method for proteomic and metabolomic data transformation (Efstathiou et al, 2017)"),
          p("",tags$b("Normalization:"),""),
          p("Median is the the default method for proteomic and metabolomic data normalization")
          )
    ),
    fluidRow(
      column(4,
             box(width = 14, title = span(tagList(icon("sliders-h"), "  ", "Setting")), status = "warning",
                 uiOutput(ns("selectSamplesUI")),
                 plotOutput(ns("completenessUI")),
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
  
  #---- sample list----
  output$selectSamplesUI <- renderUI({
    
    sampleList <- colnames(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]])
    pickerInput(
      inputId  = session$ns("selectSamples"),
      label    = .addBSpopify(label = 'Sample list:', content = "Select samples to include in further analyses"),
      choices  = sampleList,
      options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
      multiple = TRUE,
      selected = sampleList)
  })
  
  
  #---- Completeness----
  output$completenessUI <- renderPlot({
    
      plotExpDesignCompleteness(object = session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], 
                                sampleList = input$selectSamples)
  })
  
  #---- adapted parameters for each omics type ----
  output$paramUI <- renderUI({
    
    # setting for RNAseq
    paramRNAseq.list <- list(
      
      fluidRow(
        column(12, h5(.addBSpopify(label = 'Low count Filtering (CPM):', content = "Genes with low counts will be removed based on count per million (cpm), accounting for the library size"))),
        column(12,
               selectInput(inputId  = session$ns("Filter_Strategy"),
                           label    = .addBSpopify(label = 'Strategy', content = "Choose the strategy to filter genes based on count per million (cpm). keep gene if the NbOfsample_over_cpm >= Strategy."),
                           choices  = c("NbConditions" = "NbConditions",  "NbReplicates" = "NbReplicates"),
                           selected = "NbReplicates")
        ),
        column(12,
               numericInput(inputId = session$ns("FilterSeuil"),
                            label=.addBSpopify(label = 'CPM cut-off:', content = "Choose the cpm cut-off"),
                            value=1, min = 1, max=10, step = 1 )
        )
      ),
      fluidRow(
        column(12,
               h5("Gene counts normalization:"),
               selectInput(inputId  = session$ns("selectNormMethod"),
                           label    = .addBSpopify(label = 'method', content = "Normalization method"),
                           choices  =  list("TMM (edgeR)" = "TMM"),
                           selected = "TMM"))),
      fluidRow( column(12, actionButton(session$ns("run"),"Run")))
    )
    
    #---- setting for meta or RNAseq----
    paramProtMeta.list <- list(
      
      radioButtons(
        inputId  = session$ns("dataTransform"),
        label    = .addBSpopify(label = 'Data transformation:', content = "Choose log2 transformation"),
        choices  = c("log2" = "log2", "none" = "none"), 
        selected = "none"),
      hr(),
      
      radioButtons(inputId = session$ns("selectProtMetNormMethod"),
                   label    = .addBSpopify(label = 'Normalization method:', content = "Choose median normalization"),
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
      tabPanel("Distribution (boxplot)",
               tags$br(),tags$br(),
               tags$i("Expect aligned medians after running the pre-processing step. Samples with shifted median may be outliers."),
               tags$br(),tags$hr(),tags$br(),
               uiOutput(session$ns("boxplotUI"))
      ),
      tabPanel("Distribution (density)",
               tags$br(),tags$br(),
               tags$i("Expect a gaussian density distribution after running the pre-processing step. 
                      If a second peak is observed at the beginning of the curve, 
                      this could indicate that the low, uninformative
                      values, have not been correctly filtered. You may increase the filtering threshold."),
               tags$br(),tags$hr(),tags$br(),
               uiOutput(session$ns("CountDistUI"))
      ), 
      tabPanel("Principal component analysis",
               tags$br(),tags$br(),
               tags$i("You may have a look to the percentage of variability associated to each biological factor. 
               To do this, you can change the colour of samples according to experimental factor's modalities (Levels radio buttons). 
               It will help you to interpret the PCA axes. You may identify outliers samples that drives the variability. 
                Biological replicates have to group together with superposed ellipse. If not, it may indicate batch effect."),
               tags$br(),tags$hr(),tags$br(),
               uiOutput(session$ns("PCAcoordUI"))
      )
    )
    
    tabPanel.list <- list()
    
    switch(getOmicsTypes(SE.data),
           "RNAseq" = {
             tabPanel.list <- c(tabPanel.default.list,list(
                                tabPanel("Library size",
                                      tags$br(),tags$br(),
                                      tags$i("Expect equal library size after running the 
                                             pre-processing step"),
                                       tags$br(),tags$hr(),tags$br(),
                                       uiOutput(session$ns("LibSizeUI")))))
                                
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
    completeCheckRes <- checkExpDesignCompleteness(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], input$selectSamples)
    
    if(isTRUE(completeCheckRes[["error"]])){
      showModal(modalDialog(title = "Error message", completeCheckRes[["messages"]]))
    }
    validate({ need(!isTRUE(completeCheckRes[["error"]]), message=completeCheckRes[["messages"]]) })
    
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
    rea.values$datasetProcess <- rea.values$datasetProcess[rea.values$datasetProcess != dataset]
    rea.values$datasetDiff    <- rea.values$datasetDiff[rea.values$datasetDiff != dataset]
    for (database in names(rea.values$datasetDiffAnnot)){
      rea.values$datasetDiffAnnot[[database]] <- 
        rea.values$datasetDiffAnnot[[database]][rea.values$datasetDiffAnnot[[database]] != dataset]
    }
    rea.values$datasetCoEx <- rea.values$datasetCoEx[rea.values$datasetCoEx != dataset]
    for(database in names(rea.values$datasetCoExAnnot)){
      rea.values$datasetCoExAnnot[[database]] <- 
        rea.values$datasetCoExAnnot[[database]][rea.values$datasetCoExAnnot[[database]] != dataset]
    }
    
    # remove integration analysis
    session$userData$FlomicsMultiAssay <- resetFlomicsMultiAssay(session$userData$FlomicsMultiAssay, results = c("IntegrationAnalysis"))
    
    message("# 3  => Data processing: ", dataset)
   
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

