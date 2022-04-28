
##########################################
# Part3 : Data Exploratory
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
      uiOutput(ns("configUI")),
      uiOutput(ns("tabPanelUI"))
      # box(width = 4, title = "Setting", status = "warning", solidHeader = TRUE,
      #     uiOutput(ns("selectSamplesUI")),
      #     uiOutput(ns("completenessUI")),
      #     uiOutput(ns("paramUI"))
      #     ),

    )
  )
}



QCNormalizationTab <- function(input, output, session, dataset, rea.values){

  output$configUI <- renderUI({

    validate(
      need(rea.values$model != FALSE, "Please load data")
    )
    # if(rea.values$loadData == FALSE) return()

    box(width = 4, title = "Setting", status = "warning", solidHeader = TRUE,
        uiOutput(session$ns("selectSamplesUI")),
        uiOutput(session$ns("completenessUI")),
        uiOutput(session$ns("paramUI"))
        )
  })

  #### sample list  ####
  output$selectSamplesUI <- renderUI( {

    sampleList <- colnames(FlomicsMultiAssay@ExperimentList[[dataset]])

    # sample list :
    pickerInput(
      inputId  = session$ns("selectSamples"),
      label    = "Sample list :",
      choices  = sampleList,
      options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
      multiple = TRUE,
      selected = sampleList)
  })

  #### dataset completeness ####
  output$completenessUI <- renderUI({

    SE <- FlomicsMultiAssay@ExperimentList[[dataset]]
    completeCheckRes <<- CheckExpDesignCompleteness(session$userData$Design, input$selectSamples)

    if(!is.null(completeCheckRes[["error"]])){
      showModal(modalDialog(title = "Error message", completeCheckRes[["error"]]))

    }
    # continue only if message is true or warning
    validate({ need(is.null(completeCheckRes[["error"]]) ,message="ok") })

    list(
         # plot of count per condition
         renderPlot( completeCheckRes[["plot"]])
    )
  })

  output$paramUI <- renderUI({

    paramRNAseq.list <- list(

         radioGroupButtons(inputId = session$ns("Filter_Strategy"), direction = "horizontal",
                           label = "Low count Filtering stategy (CPM) :",
                           choices = c("NbConditions" = "NbConditions",  "NbReplicates" = "NbReplicates"),
                           justified = FALSE, selected = "NbConditions"),

         numericInput(inputId = session$ns("FilterSeuil"),
                      label="CPM cutoff :",
                      value=5, min = 1, max=50, step = 1 ),

         verbatimTextOutput(session$ns("FilterResults")),

         hr(),
         selectInput(inputId  = session$ns("selectNormMethod"),
                     label    = "Normalization (method) :",
                     choices  =  list("TMM (edgeR)" = "TMM"),
                     selected = "TMM"),
         hr(),
         actionButton(session$ns("Update"),"Update")
      )

    paramProtMeta.list <- list(
        radioButtons(
          inputId  =session$ns("dataImputation"),
          label = "does the data need to be impute (no implemented yet) ?",
          choices=c("yes"="yes","no"="no"),
          selected="no"),
        hr(),

        radioButtons(
          inputId  =session$ns("dataTransform"),
          label = "Does the data need to be transformed ?",
          choices=c("log1p"="log1p","log2"="log2","log10"="log10","squareroot"="squareroot","no"="none"),
          selected="none"),
        hr(),

        actionButton(session$ns("Update"),"Update")
    )

    switch (FlomicsMultiAssay@ExperimentList[[dataset]]@metadata$omicType,
            "RNAseq"       = { paramRNAseq.list },
            "proteomics"   = { paramProtMeta.list },
            "metabolomics" = { paramProtMeta.list }
    )
  })

  # processing data with default param
  # RNAseq : normalization
  # prot/meta : transformation
  param.list <- switch( FlomicsMultiAssay@ExperimentList[[dataset]]@metadata$omicType,
              "RNAseq"       = { list(Filter_Strategy = "NbConditions", CPM_Cutoff = 5, NormMethod = "TMM") },
              "proteomics"   = { list(transform_method = "none") },
              "metabolomics" = { list(transform_method = "none") } )

  print("# 5- Compute PCA for raw counts")
  FlomicsMultiAssay@ExperimentList[[dataset]] <<-  RFLOMICS::RunPCA(FlomicsMultiAssay@ExperimentList[[dataset]])
  FlomicsMultiAssay <<- process_data(FlomicsMultiAssay = FlomicsMultiAssay, dataset = dataset,
                                     samples = colnames(FlomicsMultiAssay@ExperimentList[[dataset]]), param.list = param.list)


  ## Exploratory of Biological and Technical variability
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

        tabPanel("Principal component analysis (1/2)",
                 tags$br(),
                 fluidRow(
                   tags$br(),
                   tags$br(),
                   column(width = 12,  plotOutput(session$ns("raw.PCAcoord"))),
                   column(width = 12,  plotOutput(session$ns("norm.PCAcoord")))
                 ),
                 fluidRow(
                   column(width = 12, RadioButtonsConditionUI(session$ns("rawData"))),
                   column(width = 12,
                          tags$br(),
                          UpdateRadioButtonsUI(session$ns("rawData")),
                          tags$br(),
                          actionButton(session$ns("screenshotPCA_QC"),"Screenshot")))
        ),
        tabPanel("Principal component analysis (2/2)",
                 plotOutput(session$ns("QCdesignPCA"))
                 )
    )

    tabPanel.list <- list()
     switch(FlomicsMultiAssay@ExperimentList[[dataset]]@metadata$omicType,
            "RNAseq" = {
              tabPanel.list <- c(
                list(
                  tabPanel("Library size",
                         # library size plot
                         #column(4, plotOutput(ns("LibSize"),   height = "400%")),
                         column(12, plotOutput(session$ns("raw.LibSize"))),
                         column(12, plotOutput(session$ns("norm.LibSize")))),

                  tabPanel("Data summary",
                         column(12, uiOutput(session$ns("RNAseqSummaryUI"))) ),

                  tabPanel("Distribution (boxplot)",
                         tags$br(),
                         tags$br(),
                         column(width = 12, plotOutput(session$ns("raw.boxplot"))),
                         column(width = 12, plotOutput(session$ns("norm.boxplot"))))
                  ),
                  tabPanel.default.list)
            },
            "proteomics" = {
              tabPanel.list <- c(
                list(
                  tabPanel("Data summary",
                         column(12, uiOutput(session$ns("ProtMetaSummaryUI"))) ),

                  tabPanel("Sum of XIC",
                           # library size plot
                           #column(4, plotOutput(ns("LibSize"),   height = "400%")),
                           tags$br(),
                           tags$br(),
                           column(12, plotOutput(session$ns("raw.LibSize")))),
                  tabPanel("Distribution (boxplot)",
                           tags$br(),
                           tags$br(),
                           column(width = 12, plotOutput(session$ns("norm.boxplot"))))
                ),

                tabPanel.default.list)
            },
            "metabolomics" = {
              tabPanel.list <- c(
                list(
                  tabPanel("Data summary",
                         column(12, uiOutput(session$ns("ProtMetaSummaryUI"))) ),

                  tabPanel("Sum of XIC",
                           # library size plot
                           #column(4, plotOutput(ns("LibSize"),   height = "400%")),
                           tags$br(),
                           tags$br(),
                           column(12, plotOutput(session$ns("raw.LibSize")))),
                  tabPanel("Distribution (boxplot)",
                           tags$br(),
                           tags$br(),
                           column(width = 12, plotlyOutput(session$ns("norm.boxplot"))))
                ),
                tabPanel.default.list)
            }
    )

    # Exploratory of Biological and Technical variability
    box(width = 8, title = "Exploratory of Biological and Technical variability", solidHeader = TRUE, status = "warning",
      column(12, do.call(what = tabsetPanel, args = tabPanel.list)))
  })


  ## summary tabPanel
  # => RNAseq
  output$RNAseqSummaryUI <- renderUI({

    NbGenes             <- length(names(FlomicsMultiAssay@ExperimentList[[dataset]]))
    NbGene.rowSums.zero <- length(FlomicsMultiAssay@ExperimentList[[dataset]]@metadata$rowSums.zero)
    NbGene.low.counts   <- length(FlomicsMultiAssay@ExperimentList[[paste0(dataset, ".filtred")]]@metadata$FilteredFeatures)

    list(
    strong("Number of genes:"),
    renderPrint(NbGenes),

    strong("Number of gene 0 expression:"),
    renderPrint(NbGene.rowSums.zero),

    strong("Number of gene with low counts:"),
    renderPrint(NbGene.low.counts)
    )
  })

  # => prot or meta
  output$ProtMetaSummaryUI <- renderUI({

    NbProt <- length(names(FlomicsMultiAssay@ExperimentList[[dataset]]))
    NbProtWoutNA <- dim(na.omit(assay(FlomicsMultiAssay@ExperimentList[[dataset]])))[1]

    list(
      ## Nombre de prot
      strong("Number of Features:"),
      renderPrint(NbProt),
      # Nombre de NA
      br(),
      strong("Number of features with at least 1 NA:"),
      renderPrint(NbProt - NbProtWoutNA)
    )
  })


  #### library size plot only for RNAseq data####
  # => raw data
  output$raw.LibSize <- renderPlot({

    Library_size_barplot.plot(FlomicsMultiAssay@ExperimentList[[dataset]])
  })

  # => normalized data
  output$norm.LibSize <- renderPlot({

    Library_size_barplot.plot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
  })

  ## Try plotly for metabolomics data
  if(FlomicsMultiAssay@ExperimentList[[dataset]]@metadata$omicType=="metabolomics"){
  #### boxplot only for RNAseq data ####
  # => raw
    output$raw.boxplot <- renderPlotly({
      abundanceBoxplot(FlomicsMultiAssay@ExperimentList[[dataset]])
    })
  }else{
    output$raw.boxplot <- renderPlot({
      abundanceBoxplot(FlomicsMultiAssay@ExperimentList[[dataset]])
    })
  }

  if(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType=="metabolomics"){
  # => norm
    output$norm.boxplot <- renderPlotly({
      abundanceBoxplot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
    })
  }else{
    output$norm.boxplot <- renderPlot({
      abundanceBoxplot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
    })
  }

  #### abundance distribution ####
  # => raw for all data
  output$raw.CountDist <- renderPlot(height = 300, {

    Data_Distribution_Density.plot(FlomicsMultiAssay@ExperimentList[[dataset]])
  })

  # => display only for normalized RNAseq data
  if (FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType == "RNAseq"){

      output$norm.CountDist <- renderPlot(height = 300, {
        Data_Distribution_Density.plot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
      })
  }

  #### PCA analysis ####
  # select PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "rawData")

  # select factors for color PCA plot
  callModule(RadioButtonsCondition, "rawData")




  #### PCA plot  ####
  # => raw
  output$raw.PCAcoord <- renderPlot({
    #FlomicsMultiAssay@ExperimentList[[dataset]] <<-  RunPCA(FlomicsMultiAssay@ExperimentList[[dataset]])
    PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    condGroup <- input$`rawData-condColorSelect`[1]

    RFLOMICS::plotPCA(FlomicsMultiAssay@ExperimentList[[dataset]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)

  })


  # => display only for normalized RNAseq data
  if (FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType == "RNAseq"){
        output$norm.PCAcoord <- renderPlot({
          #FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<-  RunPCA(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
          PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
          PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
          condGroup <- input$`rawData-condColorSelect`[1]

          RFLOMICS::plotPCA(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]], PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
        })
  }


  # # save current PCA plot with fixed axis & color
  # ## screenShot
  # observeEvent(input$screenshotPCA_QC, {
  #
  #   PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
  #   PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
  #   condGroup <- input$`rawData-condColorSelect`[1]
  #
  #   # file.copy(file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
  #   #           file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
  #   #           overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
  # })

  #### PCA analysis QCdesign ####
  output$QCdesignPCA <- renderPlot({

       # if(is.null(FlomicsMultiAssay@ExperimentList[[dataset]][["PCAlist"]][["raw"]])){
       #   FlomicsMultiAssay@ExperimentList[[dataset]] <<-  RunPCA(FlomicsMultiAssay@ExperimentList[[dataset]])
       # }

       mvQCdesign(FlomicsMultiAssay,data=dataset,PCA="raw", axis=5,
               pngFile=file.path(tempdir(), paste0(dataset,"_PCAdesignCoordRaw.png")))
    })

  # #### PCA analysis QCdata ####
  # output$QCdata <- renderPlot({
  #
  #   mvQCdata(FlomicsMultiAssay,data=dataset,PCA="raw",axis=5,
  #            pngFile=file.path(tempdir(), paste0(dataset,"_PCAmetaCorrRaw.png")))
  # })




  observeEvent(input$Update, {

    rea.values[[dataset]]$diffAnal  <- FALSE
    rea.values[[dataset]]$coExpAnal <- FALSE
    rea.values[[dataset]]$diffAnnot <- FALSE
    rea.values[[dataset]]$diffValid <- FALSE

    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal <<- list()
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$CoExpAnal   <<- list()
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpEnrichAnal  <<- list()
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$CoExpEnrichAnal  <<- list()


    print("# 6.bis => Update data processing...")
    #FlomicsMultiAssay <<- resetFlomicsMultiAssay(object=FlomicsMultiAssay, results=c("DiffExpAnal", "CoExpAnal", "EnrichAnal"))

    #rea.values$dataAnalysis <- TRUE

    # update processing data
    # RNAseq : normalization
    # prot/meta : transformation
    param.list <- list()
    switch( FlomicsMultiAssay@ExperimentList[[dataset]]@metadata$omicType,
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

    FlomicsMultiAssay <<- process_data(FlomicsMultiAssay = FlomicsMultiAssay, dataset = dataset, samples = input$selectSamples, param.list = param.list)

    output$FilterResults <- renderPrint({

      paste0( length(FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$FilteredFeature),
              " features filtered (from ", dim(FlomicsMultiAssay[[dataset]])[1], ")")
    })


    # => normalized data
    output$norm.LibSize <- renderPlot({

      Library_size_barplot.plot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
    })

    #### abundance distribution ####
    output$norm.CountDist <- renderPlot(height = 300, {

      Data_Distribution_Density.plot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])
    })

    #### Boxplot of distribution of normalized abundance
    output$norm.boxplot <- renderPlot({
      abundanceBoxplot(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])

    })

    #### PCA plot
    output$norm.PCAcoord <- renderPlot({

      #FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<-  RunPCA(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])

      PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
      PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
      condGroup <- input$`rawData-condColorSelect`[1]

      RFLOMICS::plotPCA(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]],
              PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup)
    })



  }, ignoreInit = TRUE)


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



process_data <- function(FlomicsMultiAssay, dataset, samples , param.list = list(Filter_Strategy = "NbConditions", CPM_Cutoff = 1, NormMethod = "TMM", transform_method = "none")){
  print("# 6- Data processing...")
  SE <- FlomicsMultiAssay@ExperimentList[[dataset]]

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

             #### Run PCA for filtred & normalized data ####
             print("# => Compute PCA ")
             SE.processed <- RunPCA(SE.processed)
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

  FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <- SE.processed
  return(FlomicsMultiAssay)
}


