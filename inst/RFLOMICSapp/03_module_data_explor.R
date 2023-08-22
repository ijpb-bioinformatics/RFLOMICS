
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


### server
QCNormalizationTab <- function(input, output, session, dataset, rea.values){
  
  #we must select contrast before
  validate(
    need(rea.values$analysis != FALSE, "Please run data processing")
  )
  
  
  local.rea.values <- reactiveValues(dataset.processed.SE = NULL, 
                                     dataset.raw.SE = NULL,
                                     compCheck = TRUE,
                                     message = NULL)
  
  
  ### sample list  ###
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
  
  # ### check completeness ###
  # completeCheckRes <- CheckExpDesignCompleteness(object = session$userData$FlomicsMultiAssay,
  #                                                sampleList = colnames(session$userData$FlomicsMultiAssay[[dataset]]))
  # # stock message in MAE
  # session$userData$FlomicsMultiAssay@metadata$completeCheck[["error"]]   <- completeCheckRes[["error"]]
  # session$userData$FlomicsMultiAssay@metadata$completeCheck[["warning"]] <- completeCheckRes[["warning"]]
  # # reactive values
  # rea.values[[dataset]]$compCheck <- TRUE
  # if(!is.null(completeCheckRes[["error"]])){ rea.values[[dataset]]$compCheck <- FALSE }
  
  
  output$completenessUI <- renderUI({
    
    completeCheckRes <- CheckExpDesignCompleteness(session$userData$FlomicsMultiAssay, input$selectSamples)
    local.rea.values$message   <- completeCheckRes[["error"]]
    
    #local.rea.values$compCheck <- TRUE
    
    if(!is.null(completeCheckRes[["error"]])){
      #local.rea.values$compCheck <- FALSE
      
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
      fluidRow( column(12, actionButton(session$ns("run"),"Run")))
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
        inputId  = session$ns("dataTransform"),
        label    = "Does the data need to be transformed?",
        choices  = c("log1p" = "log1p","squareroot" = "squareroot","none" = "none", "log2" = "log2", "log10" = "log10"), 
        selected = "none"),
      hr(),
      
      radioButtons(inputId = session$ns("selectProtMetNormMethod"),
                  label    = "method",
                  choices  =   c("median" = "median","totalSum" = "totalSum", "none" = "none"),
                  selected = "none"),
      
      
      actionButton(session$ns("run"),"Run")
    )
    
    switch (session$userData$FlomicsMultiAssay[[dataset]]@metadata$omicType,
            "RNAseq"       = { paramRNAseq.list },
            "proteomics"   = { paramProtMeta.list },
            "metabolomics" = { paramProtMeta.list }
    )
  })
  
  #### PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "factors")
  # select factors for color PCA plot
  #callModule(RadioButtonsCondition, "factors", typeFact = c("Bio", "batch"))
  
  #### PCA for metadata axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "meta")
  # select factors for color PCA plot
  #callModule(RadioButtonsCondition, "meta", typeFact = c("Meta"))
  
  #### dataset filterin summary
  output$filtSummary1UI <- renderUI({
    
    tagList(
      
      box(title = length(names(session$userData$FlomicsMultiAssay[[dataset]])), 
          width = 6, background = "maroon", 
          paste0("Initial number of ", omics.dic[[rea.values[[dataset]]$omicsType]][["variableName"]])
      ),
      box( title = length(colnames(session$userData$FlomicsMultiAssay[[dataset]])),
           width = 6, background = "light-blue", "Initial number of samples"
      )
    )
    
  })
  output$filtSummary2UI <- renderUI({
    
    if(rea.values[[dataset]]$process == FALSE) return()
    
    tagList(
      box(
        title = length(names(local.rea.values$dataset.processed.SE)), width = 6, background = "fuchsia", 
        paste0("Number of filtred ", omics.dic[[rea.values[[dataset]]$omicsType]][["variableName"]])
      ),
      box(title = length(colnames(local.rea.values$dataset.processed.SE)),
          width = 6, background = "purple", "Number of filtred samples"
      )
    )
    
  })
  
  
  #### Exploratory of Biological and Technical variability
  output$tabPanelUI <- renderUI({
    
    if(rea.values$model == FALSE) return()
    
    factors.bio   <- names(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("Bio")])
    factors.batch <- names(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("batch")])
    factors.meta  <- names(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("Meta")])
    
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
                 column(width = 6, 
                        #RadioButtonsConditionUI(session$ns("factors"))
                        radioButtons(inputId = session$ns("PCA.factor.condition"),
                                     label = 'Levels:',
                                     choices = c("groups",factors.bio, factors.batch),
                                     selected = "groups")),
                 column(width = 6, UpdateRadioButtonsUI(session$ns("factors")))
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
    
    if(any(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("Meta")))
      tabPanel.list <- c(tabPanel.list,
                         list(tabPanel("PCA for metadata",
                                       tags$br(),
                                       column(width = 12,  plotOutput(session$ns("raw.PCA.meta"))),
                                       hr(),
                                       fluidRow(
                                         column(width = 6, 
                                                #RadioButtonsConditionUI(session$ns("meta"))
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
  
  #### QC plot for raw/processed data
  # library size plot only for RNAseq data
  output$raw.LibSize <- renderPlot({
    Library_size_barplot.plot(session$userData$FlomicsMultiAssay[[dataset]], raw = TRUE)
  })
  output$norm.LibSize <- renderPlot({
    
    if(rea.values[[dataset]]$process == FALSE) return()
    Library_size_barplot.plot(session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]])
  })
  
  # value (count/intensity) distribution (boxplot/density)
  output$raw.boxplot <- renderPlot({
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[dataset]], plot = "boxplot", raw = TRUE)
  })
  output$norm.boxplot <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]], plot = "boxplot", raw = FALSE)
  })
  
  output$raw.CountDist <- renderPlot({
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[dataset]], plot = "density", raw = TRUE)
  })
  output$norm.CountDist <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    Data_Distribution_plot(session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]], plot = "density", raw = FALSE)
  })
  
  # PCA plot
  output$raw.PCAcoord <- renderPlot({
    
    PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
    PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
    #condGroup <- input$`factors-condColorSelect`[1]
    condGroup <- input$PCA.factor.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[dataset]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
    
  })
  output$norm.PCAcoord <- renderPlot({
    if(rea.values[[dataset]]$process == FALSE) return()
    
    PC1.value <- as.numeric(input$`factors-Firstaxis`[1])
    PC2.value <- as.numeric(input$`factors-Secondaxis`[1])
    #condGroup <- input$`factors-condColorSelect`[1]
    condGroup <- input$PCA.factor.condition
    
    RFLOMICS::plotPCA(local.rea.values$dataset.processed.SE, PCA = "norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
    
  })
  
  # PCA meta plot
  output$raw.PCA.meta <- renderPlot({
    
    PC1.value <- as.numeric(input$`meta-Firstaxis`[1])
    PC2.value <- as.numeric(input$`meta-Secondaxis`[1])
    #condGroup <- input$`meta-condColorSelect`[1]
    condGroup <- input$PCA.meta.condition
    
    RFLOMICS::plotPCA(session$userData$FlomicsMultiAssay[[dataset]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
  })
  
  output$norm.PCA.meta <- renderPlot({
    
    if(rea.values[[dataset]]$process == FALSE) return()
    
    PC1.value <- as.numeric(input$`meta-Firstaxis`[1])
    PC2.value <- as.numeric(input$`meta-Secondaxis`[1])
    #condGroup <- input$`meta-condColorSelect`[1]
    condGroup <- input$PCA.meta.condition
    
    RFLOMICS::plotPCA(local.rea.values$dataset.processed.SE, PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
  })
  
  #### run proprocessing : Normalisation/transformation, filtering...
  observeEvent(input$run, {
    
    if(!is.null(local.rea.values$message)){
      showModal(modalDialog(title = "Error message", local.rea.values$message))
    }
    validate({ need(is.null(local.rea.values$message), message=local.rea.values$message) })
    
    rea.values[[dataset]]$compCheck <- TRUE
    
    # re-initialize reactive values
    rea.values[[dataset]]$process   <- FALSE
    rea.values[[dataset]]$diffAnal  <- FALSE
    rea.values[[dataset]]$coExpAnal <- FALSE
    rea.values[[dataset]]$diffAnnot <- FALSE
    rea.values[[dataset]]$diffValid <- FALSE
    
    # re-initialize MAE object
    if (dataset %in% rea.values$datasetDiff){
      rea.values$datasetDiff <- rea.values$datasetDiff[-which(rea.values$datasetDiff == dataset)]
    }
    
    # new selected param
    param.list <- list()
    switch( rea.values[[dataset]]$omicsType,
            "RNAseq" = {
              param.list <- list(Filter_Strategy = input$Filter_Strategy, CPM_Cutoff = input$FilterSeuil, NormMethod=input$selectNormMethod)
            },
            "proteomics" = {
              param.list <- list(transform_method = input$dataTransform, NormMethod=input$selectProtMetNormMethod)
            },
            "metabolomics" = {
              param.list <- list(transform_method = input$dataTransform, NormMethod=input$selectProtMetNormMethod)
            }
    )
    
    print(paste0("# 3  => Data processing : ", dataset))
    processed.SE <- process_data(SE = session$userData$FlomicsMultiAssay[[dataset]], dataset = dataset, 
                                 samples = input$selectSamples, param.list = param.list)
    
    local.rea.values$dataset.processed.SE <- processed.SE
    
    ## add new SE with processed data
    SE.name <- paste0(dataset,".filtred")
    
    # remove SE processed if existe
    if (SE.name %in% names(session$userData$FlomicsMultiAssay)){
      session$userData$FlomicsMultiAssay <- session$userData$FlomicsMultiAssay[,, -which(names(session$userData$FlomicsMultiAssay) == SE.name)]
    }
    
    session$userData$FlomicsMultiAssay <- eval(parse(text = paste0('c( session$userData$FlomicsMultiAssay ,', SE.name, ' = processed.SE )')))
    
    rea.values[[dataset]]$process   <- TRUE
    
  }, ignoreInit = TRUE)
  
  return(input)
  
}


############## functions ###############
###
process_data <- function(SE, dataset, samples , param.list = list(Filter_Strategy = "NbConditions", 
                                                                  CPM_Cutoff = 1, NormMethod = "TMM", transform_method = "none")){
  
  print("#    => select samples")
  SE.new <- SE[, SE$primary %in% samples]
  
  SE.new@metadata$Groups <- dplyr::filter(SE@metadata$Groups, samples %in% SE.new$primary)
  
  switch(SE.new@metadata$omicType,
         "RNAseq" = {
           #### Filter low abundance ####
           print("#    => Low counts Filtering...")
           SE.processed <- FilterLowAbundance(SE.new, param.list[["Filter_Strategy"]], param.list[["CPM_Cutoff"]])
           
           #### Run Normalisation ####
           print("#    => Counts normalization...")
           SE.processed <- RunNormalization(SE.processed, NormMethod = param.list[["NormMethod"]])
         },
         "proteomics" = {
           print("#    => transformation data...")
           SE.processed <- TransformData(SE.new, transform_method = param.list[["transform_method"]])
           
           if(param.list[["NormMethod"]] != "none"){
             print("#    => Run normalization...")
             SE.processed <- RunNormalization(SE.processed, NormMethod = param.list[["NormMethod"]])
           }else{SE.processed@metadata[["Normalization"]]$methode <- "none" }
           
         },
         "metabolomics" = {
           print("#    => transformation data...")
           SE.processed <- TransformData(SE.new, transform_method = param.list[["transform_method"]])
           
           if(param.list[["NormMethod"]] != "none"){
             print("#    => Run normalization...")
             SE.processed <- RunNormalization(SE.processed, NormMethod = param.list[["NormMethod"]])
           }else{SE.processed@metadata[["Normalization"]]$methode <- "none" }
         }
  )
  
  #### Run PCA for filtred & normalized data ####
  print("#    => Compute PCA ")
  SE.processed <- RunPCA(SE.processed) # 19/04/23 : nothing to change here, already transformed and replaced, normmethod available. 
  
  return(SE.processed)
}


