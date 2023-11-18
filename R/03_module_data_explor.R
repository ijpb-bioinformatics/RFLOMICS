
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
          width = 12, status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          
          p("For each diagnostic plot, both raw and processed (filtered, normalized, ...) data are displayed with expertised default parameters."),
          p("- You may first have a look at the default processed plots to eventually identify outliers sample"),
          p("- You may check onto the PCA factorial map to control that samples are groupped with the replicat of the same condition."),
          p("- You may play with filtered parameters to see if it improves the grouping"),
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
    completeCheckRes <- CheckExpDesignCompleteness(session$userData$FlomicsMultiAssay, paste0(dataset, ".raw"), input$selectSamples)

    if(isTRUE(completeCheckRes[["error"]])){

      local.rea.values$message   <- completeCheckRes[["summary"]][1,3]
    }
    
    list(
      # plot of count per condition
      renderPlot(completeCheckRes[["plot"]])
    )
  })
  
  #---- adapted parameters for each omics type ----
  output$paramUI <- renderUI({
    
    # setting for RNAseq
    paramRNAseq.list <- list(
      
      fluidRow(
        column(12, h5("Low count Filtering (CPM):")),
        column(8,
               selectInput(inputId  = session$ns("Filter_Strategy"),
                           label    = "Stategy",
                           choices  = c("NbConditions" = "NbConditions",  "NbReplicates" = "NbReplicates"),
                           selected = "NbConditions")
        ),
        column(4,
               numericInput(inputId = session$ns("FilterSeuil"),
                            label="Cutoff:",
                            value=5, min = 1, max=10, step = 1 )
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
    
    factors.bio   <- bioFactors(MAE.data)
    factors.batch <- batchFactors(MAE.data)
    factors.meta  <- metaFactors(MAE.data)
    
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
                 column(width = 6, 
                        radioButtons(inputId = session$ns("PCA.factor.condition"),
                                     label = 'Levels:',
                                     choices = c("groups", factors.bio, factors.batch),
                                     selected = "groups")),
                 column(width = 6, UpdateRadioButtonsUI(session$ns("factors")))
               ),
               hr(),
               column(width = 12, plotOutput(session$ns("norm.PCAcoord")))
      )
    )
    
    tabPanel.list <- list()
    switch(getOmicsTypes(SE.data),
           "RNAseq" = {
             tabPanel.list <- c(
               list(
                 tabPanel("Library size",
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
    
    if (any(getFactorTypes(MAE.data) %in% c("Meta")))
      tabPanel.list <- c(tabPanel.list,
                         list(tabPanel("PCA for metadata",
                                       tags$br(),
                                       column(width = 12,  plotOutput(session$ns("raw.PCA.meta"))),
                                       hr(),
                                       fluidRow(
                                         column(width = 6, 
                                                radioButtons(inputId = session$ns("PCA.meta.condition"),
                                                             label = 'Levels:',
                                                             choices = factors.meta,
                                                             selected = factors.meta[1])),
                                         column(width = 6, UpdateRadioButtonsUI(session$ns("meta")))
                                       ),
                                       hr(),
                                       column(width = 12, plotOutput(session$ns("norm.PCA.meta")))
                         )))
    
    # Exploratory of Biological and Technical variability
    box(width = 14, title = "Exploratory of Biological and Technical variability", 
        solidHeader = TRUE, status = "warning",
        do.call(what = tabsetPanel, args = tabPanel.list))
  })
  
  #---- QC plot for raw/processed data----
  # library size plot only for RNAseq data
  output$raw.LibSize <- renderPlot({
    Library_size_barplot.plot(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], raw = TRUE)
  })
  
  output$norm.LibSize <- renderPlot({
    
    if(rea.values[[dataset]]$process == FALSE) return()
    Library_size_barplot.plot(session$userData$FlomicsMultiAssay[[dataset]])
  })
  
  # value (count/intensity) distribution (boxplot/density)
  output$raw.boxplot <- renderPlot({
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], plot = "boxplot", raw = TRUE)
  })
  output$norm.boxplot <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[dataset]], plot = "boxplot", raw = FALSE)
  })
  
  output$raw.CountDist <- renderPlot({
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], plot = "density", raw = TRUE)
  })
  output$norm.CountDist <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[dataset]], plot = "density", raw = FALSE)
  })
  
  # PCA plot
  output$raw.PCAcoord <- renderPlot({
    
    PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
    PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
    condGroup <- input$PCA.factor.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
    
  })
  output$norm.PCAcoord <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
    PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
    condGroup <- input$PCA.factor.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[dataset]], PCA = "norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
    
  })
  
  # PCA meta plot
  output$raw.PCA.meta <- renderPlot({
    
    PC1.value <- as.numeric(input$`meta-Firstaxis`[1])
    PC2.value <- as.numeric(input$`meta-Secondaxis`[1])
    condGroup <- input$PCA.meta.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
  })
  
  output$norm.PCA.meta <- renderPlot({
    
    if(rea.values[[dataset]]$process == FALSE) return()
    
    PC1.value <- as.numeric(input$`meta-Firstaxis`[1])
    PC2.value <- as.numeric(input$`meta-Secondaxis`[1])
    condGroup <- input$PCA.meta.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[dataset]], PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
  })
  
  #---- run preprocessing - Normalization/transformation, filtering...----
  observeEvent(input$run, {
    
    # check completeness for curent dataset
    if(!is.null(local.rea.values$message)){
      showModal(modalDialog(title = "Error message", local.rea.values$message))
    }
    validate({ need(is.null(local.rea.values$message), message=local.rea.values$message) })
    
    
    # get parameters 
    param.list <- list()
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
    
    print(paste0("# 3  => Data processing: ", dataset))
    processed.SE <- runDataProcessing(session$userData$FlomicsMultiAssay[[paste0(dataset, ".raw")]], samples = input$selectSamples,  
                                     lowCountFiltering=list(strategy=param.list[["Filter_Strategy"]], CPM_Cutoff=param.list[["CPM_Cutoff"]]), 
                                     normalisation=list(method=param.list[["NormMethod"]]), transformation=list(method=param.list[["transform_method"]]))
    
    
    # remove SE processed if exist
    if (dataset %in% names(session$userData$FlomicsMultiAssay)){
      session$userData$FlomicsMultiAssay <- session$userData$FlomicsMultiAssay[,, -which(names(session$userData$FlomicsMultiAssay) == dataset)]
    }
    
    # add new SE with processed data
    session$userData$FlomicsMultiAssay <- eval(parse(text = paste0('c( session$userData$FlomicsMultiAssay ,', dataset, ' = processed.SE )')))
    
    rea.values[[dataset]]$process <- TRUE
    
  }, ignoreInit = TRUE)
  
  return(input)
  
}


############## functions ###############

# ----- check run norm execution ------
check_run_process_execution <- function(object.MAE, dataset, param.list = NULL){

  if(!dataset %in% names(object.MAE))  return(TRUE)

  SE <- object.MAE[[dataset]]
  switch( getOmicsTypes(object.MAE[[dataset]]),
          "RNAseq" = {
            # filtering setting
            if(param.list$Filter_Strategy != getFilterSetting(SE)$filterStrategy) return(TRUE)
            if(param.list$CPM_Cutoff      != getFilterSetting(SE)$cpmCutoff)      return(TRUE)
            
            # normalisation setting
            if(param.list$NormMethod      != getNormSetting(SE)$method) return(TRUE)
          },
          "proteomics" = {

            if(param.list$transform_method != getTransSetting(SE)$method) return(TRUE)
            if(param.list$NormMethod       != getNormSetting(SE)$method) return(TRUE)
          },
          "metabolomics" = {

            if(param.list$transform_method != getTransSetting(SE)$method) return(TRUE)
            if(param.list$NormMethod       != getNormSetting(SE)$method) return(TRUE)
          }
  )

  return(FALSE)
}

