
#' @importFrom UpSetR upset
#' @importFrom DT formatSignif datatable renderDataTable styleInterval 
#' formatStyle styleEqual
#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput 
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom dplyr filter select arrange
#' @importFrom purrr reduce
#' @importFrom magrittr "%>%"

DiffExpAnalysisUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    fluidRow(
      uiOutput(ns("instruction"))),
    
    ### parametres for Diff Analysis
    fluidRow(
      column(3,
             uiOutput(ns("DiffParamUI"))),
      column(9,
             uiOutput(ns("ContrastsResults")),
             uiOutput(ns("validateUI")))),
    
    fluidRow(
      uiOutput(ns("ResultsMerge")))
  )
}



DiffExpAnalysis <- function(input, output, session, dataset, rea.values){
  
  local.rea.values <- reactiveValues(p.adj.cutoff = 0.05, abs.logFC.cutoff = 0, selectedContrasts = NULL)
  
  # list of tools for diff analysis
  MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")
  
  method <- switch (rea.values[[dataset]]$omicsType,
                    "RNAseq"       = MethodList[1],
                    "proteomics"   = MethodList[2],
                    "metabolomics" = MethodList[2])
  
  output$instruction <- renderUI({
    box(title = span(tagList(icon("cogs"), "  ",  a(names(method), href="https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf"), tags$small("(Scroll down for instructions)")  )),
        solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
        p("Differential expression analysis is performed for each contrast. 
          There are just two options to set (the adjusted-pvalue cut-off and the |logFC| cut-off).
          The results will appear in blocks with the contrast's name and statistics (one per contrast), each block offering a tab panel with several outputs:"),
        p("- The graph of Pvalue's distribution: Distribution of pvalue's which has to be check to validate results. The most desirable shape is a pick of p-values at 0 following by a uniform distribution. ", a("(some help to identify the good shapes)", href="/www/Pvalue_distrib.pdf"),""),
        p("- The MA plot which gives the logFC across the mean of the expression/abundance"),
        p("- The Volcano plot: implemented in the EnhancedVolcano R-package (Blighe K, Rana S, Lewis M (2022). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling R package version 1.12.0.))"),
        p("- A dataframe with the statistical results of the differential analysis per DE entities"),
        p("- A Heatmap plot: implemented in the ComplexHeatmap package (Gu Z (2022). Complex Heatmap Visualization. iMeta. doi:10.1002/imt2.43.)"),
        p("- A PCA on DE entities"),
        p("- Boxplot of DE : boxplot showing the expression/abundance profile of a selected DE entity across experimental factors"),
    )
  })
  
  # diff param
  output$DiffParamUI <- renderUI({
    
    #we must run process before
    validate(
      need(rea.values[[dataset]]$process != FALSE, "Please run data processing")
    )
    
    validate(
      need(!is.null(rea.values$Contrasts.Sel), "Please run data processing"))
    
    #design must be complete
    validate(
      need(rea.values[[dataset]]$compCheck != FALSE, 
           metadata(session$userData$FlomicsMultiAssay)$completeCheck[["error"]])
    )
    
    local.rea.values$selectedContrasts <- 
      getSelectedContrasts(session$userData$FlomicsMultiAssay[[dataset]])
      # generateExpressionContrast(session$userData$FlomicsMultiAssay[[dataset]]) %>% 
      # reduce(rbind) %>% 
      # filter(contrast %in% rea.values$Contrasts.Sel$contrast)
    
    validate(
      need(nrow(local.rea.values$selectedContrasts) != 0, 
           message = "no contrast matches the sample selection")
    )
    
    
    ## getcontrast
    box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
        
        fluidRow(column(12,
                        
                        ## list of contrasts to test
                        pickerInput(
                          inputId  = session$ns("contrastList"),
                          label    = .addBSpopify(label = 'Selected contrasts:', content = "Contrasts to run diff analysis. If you want to test all contrasts, select 'All'"),
                          choices  = local.rea.values$selectedContrasts$contrastName),
                        
                        # method for Diff analysis
                        selectInput(inputId  = session$ns("AnaDiffMethod"),
                                    label = .addBSpopify(label = 'Method:', content = "Differential analysis  method"),
                                    choices  = method,
                                    selected = method),
                        
                        numericInput(inputId = session$ns("p.adj.cutoff"), 
                                     label=.addBSpopify(label = 'Adjusted pvalue cutoff:', content = "The adjusted p-value cut-off"),
                                     value=local.rea.values$p.adj.cutoff, min=0, max=1, 0.01),
                        numericInput(inputId = session$ns("abs.logFC.cutoff"),
                                     label=.addBSpopify(label = '|logFC| cutoff:', content = "the absolute log FC cut-off"),
                                     value=local.rea.values$abs.logFC.cutoff, min=0, max=100, 0.01),
                        
                        # use of cluster. need setting step
                        materialSwitch(inputId = session$ns("clustermq"),
                                       label=.addBSpopify(label = 'use remote Cluster:', content = "send calculation to the cluster"),
                                       value = FALSE, status = "success"),
                        
                        actionButton(session$ns("runAnaDiff"),"Run"))
        ))
    
  })
  
  # # filter param
  output$validateUI <- renderUI({
    
    if (rea.values[[dataset]]$diffAnal == FALSE) return()
    if (is.null(rea.values[[dataset]]$DiffValidContrast) || dim(rea.values[[dataset]]$DiffValidContrast)[1] == 0) return()
    
    
    fluidRow(
      column(width = 9),
      column(width = 3, actionButton(session$ns("validContrast"),"Validate")))
  })
  
  # Run the differential analysis for each contrast set
  # Filter
  #   -> return a dynamic user interface with a collapsible box for each contrast
  #         - Pvalue graph
  #         - MAplot
  #         - Table of the DE genes
  #   -> combine data : union or intersection
  
  #### PCA axis for plot
  observe({
    lapply(seq_len(length(rea.values$Contrasts.Sel$contrast)), function(i) {
      
      vect     <- unlist(rea.values$Contrasts.Sel[i,])
      
      # update/adapt PCA axis
      callModule(UpdateRadioButtons, paste0(vect["contrastName"],"-diff"))
    })
  })
  
  ##================================ RUN =======================================##
  
  ### run diff
  #######################
  observeEvent(input$runAnaDiff, {
    
    param.list <- list(method             = input$AnaDiffMethod,
                       clustermq          = input$clustermq,
                       p.adj.method  = "BH",
                       p.adj.cutoff  = input$p.adj.cutoff, 
                       abs.logFC.cutoff   = input$abs.logFC.cutoff)
    
    if(check_run_diff_execution(session$userData$FlomicsMultiAssay[[dataset]], param.list) == FALSE) return()
    
    rea.values[[dataset]]$diffAnal   <- FALSE
    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    
    session$userData$FlomicsMultiAssay <-
      resetRflomicsMAE(session$userData$FlomicsMultiAssay,
                       datasetNames = dataset,
                       singleAnalyses = c("DiffExpAnal",
                                    "DiffExpEnrichAnal",
                                    "CoExpAnal",
                                    "CoExpEnrichAnal"),
                       multiAnalyses = c("IntegrationAnalysis"))
    
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$DiffExpEnrichAnal <- list()
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$CoExpAnal         <- list()
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$CoExpEnrichAnal   <- list()
    
    # reset reactive values
    rea.values$datasetDiff <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "DiffExpAnal")
    rea.values$datasetDiffAnnot <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "DiffExpEnrichAnal")
    rea.values$datasetCoEx <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "CoExpAnal")
    rea.values$datasetCoExAnnot <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "CoExpEnrichAnal")
    
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$DiffExpAnal[["Validcontrasts"]] <- NULL
    
    # session$userData$FlomicsMultiAssay[[dataset]] <- 
    #   setValidContrasts(session$userData$FlomicsMultiAssay[[dataset]], 
    #                     contrastList = NULL)
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    if(length(metadata(dataset.SE)$DiffExpAnal) == 0){
      
      message("# 4- Diff Analysis... ", dataset)
      message("#    => Filter Diff Analysis...")
      
      # run diff analysis with selected method
      dataset.SE <- runDiffAnalysis(
        object        = dataset.SE,
        p.adj.method  = "BH",
        method        = input$AnaDiffMethod,
        clustermq     = input$clustermq,
        p.adj.cutoff  = input$p.adj.cutoff, 
        logFC.cutoff  = input$abs.logFC.cutoff,
        cmd           = TRUE)
    }
    else{
      
      message("#    => Filter Diff Analysis...")
      ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
      dataset.SE <- filterDiffAnalysis(
        object        = dataset.SE,
        p.adj.cutoff  = input$p.adj.cutoff,
        logFC.cutoff  = input$abs.logFC.cutoff)
    }
    
    # error management
    if(isFALSE(metadata(dataset.SE)$DiffExpAnal[["results"]])){
      showModal(
        modalDialog(
          title = "Error message",
          if(!is.null(metadata(dataset.SE)$DiffExpAnal[["ErrorStats"]])){
            DT::renderDataTable(
              metadata(dataset.SE)$DiffExpAnal[["ErrorStats"]],
              rownames = FALSE)
          }
          else{
            as.character(metadata(dataset.SE)$DiffExpAnal[["Error"]])
          }
        ))
    }
    
    if(is.null(metadata(dataset.SE)$DiffExpAnal[["RawDEFres"]])){
      
      showModal(
        modalDialog(
          title = "Error message",
          if(!is.null(metadata(dataset.SE)$DiffExpAnal[["ErrorTab"]])){
            DT::renderDataTable(
              metadata(dataset.SE)$DiffExpAnal[["ErrorTab"]],
              rownames = FALSE)
          }
          else{
            as.character(metadata(dataset.SE)$DiffExpAnal[["error"]])
          }
        ))
    }
    
    session$userData$FlomicsMultiAssay[[dataset]] <- dataset.SE
    
    rea.values[[dataset]]$diffAnal <- TRUE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#
    
  }, ignoreInit = TRUE)
  
  ### validate contrasts
  #######################
  observeEvent(input$validContrast, {
    
    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    
    rea.values$datasetDiff <-
      rea.values$datasetDiff[rea.values$datasetDiff != dataset]
    
    session$userData$FlomicsMultiAssay <-
      resetRflomicsMAE(session$userData$FlomicsMultiAssay,
                       datasetNames = dataset,
                       singleAnalyses = c("DiffExpEnrichAnal",
                                          "CoExpAnal",
                                          "CoExpEnrichAnal"),
                       multiAnalyses = c("IntegrationAnalysis"))
    
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$DiffExpEnrichAnal <- NULL
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$CoExpEnrichAnal   <- NULL
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$CoExpAnal         <- NULL
    
    # metadata(session$userData$FlomicsMultiAssay[[dataset]])$DiffExpAnal[["Validcontrasts"]] <- 
    #   rea.values[[dataset]]$DiffValidContrast
    
    session$userData$FlomicsMultiAssay[[dataset]] <- 
      setValidContrasts(session$userData$FlomicsMultiAssay[[dataset]],
                        contrastList = rea.values[[dataset]]$DiffValidContrast)
    
    rea.values[[dataset]]$diffValid <- TRUE
    
    # reset reactive values
    rea.values$datasetDiff <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "DiffExpAnal")
    rea.values$datasetDiffAnnot <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "DiffExpEnrichAnal")
    rea.values$datasetCoEx <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "CoExpAnal")
    rea.values$datasetCoExAnnot <-
      getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                              analyses = "CoExpEnrichAnal")
    
  }, ignoreInit = TRUE)
  
  ##============================== DISPLAY =====================================##
  
  # display results per contrast
  output$ContrastsResults <- renderUI({
    
    if (rea.values[[dataset]]$diffAnal == FALSE ||
        is.null(metadata(session$userData$FlomicsMultiAssay[[dataset]])$DiffExpAnal[["TopDEF"]])) return()
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    # list of bio factors
    factors.bio   <- getBioFactors(dataset.SE)
    factors.batch <- getBatchFactors(dataset.SE)
    factors.meta  <- getMetaFactors(dataset.SE)
    
    list(
      lapply(seq_len(nrow(getSelectedContrasts(dataset.SE))), function(i) {
        
        vect     <- unlist(metadata(dataset.SE)$design$Contrasts.Sel[i,])
        res      <- metadata(dataset.SE)$DiffExpAnal[["RawDEFres"]][[vect["contrastName"]]]
        stats    <- metadata(dataset.SE)$DiffExpAnal[["stats"]][vect["contrastName"],]
        
        diff.plots <- plotDiffAnalysis(dataset.SE, 
                                       contrastName = vect["contrastName"])
        
        if (dim(metadata(dataset.SE)$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]])[1] == 0){
          
          
          fluidRow(
            column(10,
                   box(width=14, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "danger",
                       title = tags$h5(paste0(vect["tag"], " : ", vect["contrastName"],"  [#DE: ", stats[["All"]]," ]")),
                       
                       tabsetPanel(
                         
                         ### pvalue plot ###
                         tabPanel("Pvalue's distribution", renderPlot({ diff.plots$Pvalue.hist })),
                         
                         ### MAplot
                         tabPanel("MA plot", renderPlot({ suppressMessages(diff.plots$MA.plot) })),
                         
                         ### MAplot
                         tabPanel("Volcano plot", renderPlot({ suppressMessages(diff.plots$Volcano.plot) }, height = 600))
                       )
                   )))
        }
        else{
          
          fluidRow(
            column(10,
                   box(width=14, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "success",
                       title = tags$h5(paste0(vect["tag"], ": ", vect["contrastName"],
                                              "  [#DE: ", stats[["All"]], " ; ",
                                              "Up: ",stats[["Up"]]," (",round(stats[["Up"]]/stats[["All"]],2)*100," %)"," ; ",
                                              "Down: ", stats[["Down"]]," (",round(stats[["Down"]]/stats[["All"]],2)*100,"%)]")),
                       
                       tabsetPanel(
                         
                         ### pvalue plot ###
                         tabPanel("Pvalue's distribution", 
                                  tags$br(),
                                  tags$i("You must have a look at to the distribution of non-adjusted p-values to validate
                                         your analysis. The most desirable shape is a peak of p-values at 0 following by a uniform 
                                         distribution"),
                                  tags$br(),tags$hr(),tags$br(),
                                  renderPlot({ diff.plots$Pvalue.hist })),
                         
                         ### MAplot
                         tabPanel("MA plot", 
                                  tags$br(),
                                  tags$i(paste0("Expect a majority of ", .omicsDic(dataset.SE)$variableName," around 0.",
                                                " The red points are the ", .omicsDic(dataset.SE)$variableName,
                                                " significantly over-expressed in the left factor's modality in the contrast expression whereas ",
                                                "blue points are ", .omicsDic(dataset.SE)$variableName, 
                                                " significantly under-expressed in the rigth factor's modality in the contrast expression.
                                                Only the top 20 ", .omicsDic(dataset.SE)$variableName, " DE are labeled.")),
                                  tags$br(),tags$hr(),tags$br(),
                                  renderPlot({ diff.plots$MA.plot })),
                         
                         ### MAplot
                         tabPanel("Volcano plot", 
                                  tags$br(),
                                  tags$i(paste0(" Red points are ",.omicsDic(dataset.SE)$variableName," of interest: ",
                                                "displaying both large magnitude log-fold-changes (x axis) and high statistical significance (y axis)",
                                                "Only the top 20 ", .omicsDic(dataset.SE)$variableName, " DE are labeled.")),
                                  tags$br(), tags$hr(), tags$br(),
                                  renderPlot({ diff.plots$Volcano.plot }, height = 600)),
                         
                         ### DEF result table ###
                         tabPanel("Table",
                                  tags$br(),
                                  tags$i("Table of results of the differential expression/abundance statistical analysis:"),
                                  tags$ul(
                                    tags$li(tags$i(paste0("row names: ",omicsDic(dataset.SE)$variableName," ID"))),
                                    tags$li(tags$i("logFC: log2 fold change")),
                                    tags$li(tags$i("Abundance: mean expression/abundance for the factor's modality")),
                                    tags$li(tags$i("t: t-statistic (limma-lmFit, prot/metabo)")),
                                    tags$li(tags$i("pvalue: p-values")),
                                    tags$li(tags$i("Adj.value: adjusted p-value (BH)")),
                                    tags$li(tags$i("LR: likelihood ratio test (edgeR-glmLRT, RNAseq)")),
                                    tags$li(tags$i("B: log-odds that the prot/metabo is differentially expressed (limma-topTable)")),
                                    tags$li(tags$i("Regulation = Up (green) or Down (red) regulated"))
                                  ),
                                  tags$hr(), tags$br(),
                                  ### DEF result table ###
                                  DT::renderDataTable({
                                    resTable <- metadata(dataset.SE)$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]]
                                    resTable$Regulation <- ifelse(resTable$logFC > 0, "Up", "Down")
                                    resTable %>% DT::datatable(
                                      extensions = 'Buttons',
                                      options = list(dom = 'lfrtipB',
                                                     rownames = FALSE,
                                                     pageLength = 10,
                                                     buttons = c('csv', 'excel'),
                                                     lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All")))) %>%
                                      formatStyle('Regulation',
                                                  backgroundColor = DT::styleEqual(c("Up", "Down"), 
                                                                                   c("#C7DCA7", c("#FFC5C5"))),
                                                  fontWeight = 'bold') %>% 
                                      formatSignif(columns = seq_len(dim(resTable)[2]), digits = 3)
                                  }, server = FALSE)),
                         
                         ### Heatmap ###
                         tabPanel("Heatmap",
                                  tags$br(),
                                  tags$i(tags$p(paste0("Heatmap is performed on DE ",.omicsDic(dataset.SE)$variableName,
                                                       " expression data table which has been",
                                                       " transformed by: ",getTransSettings(dataset.SE)$method, " method", 
                                                       " and normalized by: ", getNormSettings(dataset.SE)$method ,"  method."))),
                                  tags$i(tags$p(paste0("Clustering is independently performed on samples (row) and centered ",
                                                       .omicsDic(dataset.SE)$variableName,
                                                       " (column) using euclidian distance and complete aggregation method."))),
                                  tags$i(tags$p(" You may separate the heatmap by modality of factor of interest (Levels radio buttons). You may
                                         also add annotations to the heatmap by selecting biological/batch factors to display.")),
                                  tags$br(), tags$hr(), tags$br(),
                                  renderPlot({
                                    annot_arg <- c(input[[paste0(vect["contrastName"], "-annotBio")]],
                                                   input[[paste0(vect["contrastName"], "-annotBatch")]])
                                    if (length(factors.meta) > 0) {
                                      annot_arg <- c(annot_arg, input[[paste0(vect["contrastName"], "-annotMeta")]])
                                    }
                                    
                                    plotHeatmapDesign(object     = dataset.SE, 
                                                      contrastName = vect["contrastName"], 
                                                      splitFactor  = input[[paste0(vect["contrastName"],"-heat.condColorSelect")]],
                                                      annotNames =  annot_arg)
                                  }),
                                  ## select cluster to plot
                                  column(6, radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-heat.condColorSelect")),
                                                         label = 'Levels:',
                                                         choices = c("none", factors.bio),
                                                         selected = "none", inline = TRUE)),
                                  
                                  ## select annotations to show
                                  
                                  
                                  column(6 ,checkboxGroupInput(inputId = session$ns(paste0(vect["contrastName"], "-annotBio")),
                                                               label = "Biological factors", inline = TRUE,
                                                               choices = factors.bio,
                                                               selected = factors.bio)),
                                  column(6, checkboxGroupInput(inputId = session$ns(paste0(vect["contrastName"], "-annotBatch")),
                                                               label = "Batch factors",  inline = TRUE,
                                                               choices = factors.batch,
                                                               selected = NULL)),
                                  if (length(factors.meta) > 0) {
                                    column(6,checkboxGroupInput(inputId = session$ns(paste0(vect["contrastName"], "-annotMeta")),
                                                                label = "Metadata factors", inline = TRUE,
                                                                choices = factors.meta,
                                                                selected = NULL))
                                  }
                         ),
                         
                         ### PCA ###
                         tabPanel("PCA on DE",
                                  tags$br(),
                                  tags$i("PCA plot of the DE entities"),
                                  tags$br(), tags$hr(), tags$br(),
                                  fluidRow(
                                    column(width = 12,
                                           renderPlot({
                                             
                                             newDataset.SE <- dataset.SE[, which(colnames(dataset.SE) %in% dataset.SE$samples)]
                                             newDataset.SE <-  runOmicsPCA(newDataset.SE[row.names(metadata(newdataset.SE)$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]])],ncomp = 5, raw = FALSE) 
                                             
                                             PC1.value <- as.numeric(input[[paste0(vect["contrastName"],"-diff-Firstaxis")]][1])
                                             PC2.value <- as.numeric(input[[paste0(vect["contrastName"],"-diff-Secondaxis")]][1])
                                             condGroup <- input[[paste0(vect["contrastName"],"-pca.DE.condColorSelect")]][1]
                                             
                                             plotOmicsPCA(newDataset.SE,  raw = "norm", axes = c(PC1.value, PC2.value), groupColor = condGroup)
                                           })
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 6, 
                                           radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-pca.DE.condColorSelect")),
                                                        label = 'Levels:',
                                                        choices = c("groups",factors.bio),
                                                        selected = "groups")),
                                    column(width = 6, UpdateRadioButtonsUI(session$ns(paste0(vect["contrastName"],"-diff"))))
                                    
                                  )
                         ),
                         
                         ### boxplot DE ###
                         tabPanel("boxplot DE",
                                  tags$br(),
                                  tags$i(paste0("Boxplot showing the expression/abundance profile of a selected DE ",.omicsDic(dataset.SE)$variableName),
                                         " colored by experimental factor's modality (see Levels radio buttons)."),
                                  tags$br(), tags$hr(), tags$br(),
                                  fluidRow(
                                    
                                    column(width = 3,
                                           
                                           selectizeInput(
                                             inputId = session$ns(paste0(vect["contrastName"], "-DE")), 
                                             label = paste0("Select DE ",.omicsDic(dataset.SE)$variableName,":"),
                                             choices = rownames(arrange(metadata(dataset.SE)$DiffExpAnal$TopDEF[[vect["contrastName"]]], 
                                                                        Adj.pvalue)),
                                             multiple = FALSE),
                                           radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-DEcondition")),
                                                        label = 'Levels:',
                                                        choices = c("groups",factors.bio),
                                                        selected = factors.bio[1])
                                    ),
                                    column(width = 9,
                                           renderPlot({
                                             plotBoxplotDE(object=dataset.SE, features=input[[paste0(vect["contrastName"], "-DE")]], 
                                                           groupColor=input[[paste0(vect["contrastName"], "-DEcondition")]]) })
                                    )
                                  )
                         )
                       )
                   )
            ),
            column(width = 2,
                   checkboxInput(session$ns(paste0("checkbox_", vect[["tag"]])), "OK", value = TRUE)) 
          )
        }
      })
    )
    
  })
  
  # merge results on upset plot
  output$ResultsMerge <- renderUI({
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    if (rea.values[[dataset]]$diffAnal == FALSE ||
        is.null(metadata(dataset.SE)$DiffExpAnal[["mergeDEF"]])) return()
    
    DEF_mat <- as.data.frame(metadata(dataset.SE)$DiffExpAnal[["mergeDEF"]])
    
    index <- sapply(names(DEF_mat)[-1], function(x){(input[[paste0("checkbox_",x)]])}) %>% unlist()
    
    H_selected <- names(DEF_mat)[-1][index]
    
    rea.values[[dataset]]$DiffValidContrast <- filter(metadata(dataset.SE)$design$Contrasts.Sel, tag %in% H_selected)
    
    if (length(H_selected) > 1 && dim(DEF_mat)[1] != 0){
      
      box(width=12,  status = "warning",
          
          renderPlot({upset(DEF_mat, sets = H_selected, order.by = "freq") })
      )
    }
  })
  
  return(input)
}

############## functions ###############

# ----- check run diff execution ------
check_run_diff_execution <- function(object.SE, param.list = NULL){
  
  # filtering setting
  if (length(metadata(object.SE)[["DiffExpAnal"]]) == 0 || 
      metadata(object.SE)[["DiffExpAnal"]]$results != TRUE) return(TRUE)
  
  if(param.list$method != getDiffSettings(object.SE)$method) return(TRUE)
  if(param.list$p.adj.method  != getDiffSettings(object.SE)$p.adj.method) return(TRUE)
  if(param.list$p.adj.cutoff  != getDiffSettings(object.SE)$p.adj.cutoff) return(TRUE)
  if(param.list$abs.logFC.cutoff   != getDiffSettings(object.SE)$abs.logFC.cutoff)  return(TRUE)
  
  return(FALSE)
}


