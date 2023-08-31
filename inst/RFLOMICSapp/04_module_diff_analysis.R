

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
  
  local.rea.values <- reactiveValues(dataset.SE = NULL, Adj.pvalue.cutoff = 0.05, abs.logFC.cutoff = 0, update = FALSE)
  
  # list of tools for diff analysis
  MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")
  
  method <- switch (rea.values[[dataset]]$omicsType,
                    "RNAseq"       = MethodList[1],
                    "proteomics"   = MethodList[2],
                    "metabolomics" = MethodList[2])
  
  output$instruction <- renderUI({
    box(title = span(tagList(icon("cogs"), "  ",  a(names(method), href="https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf"), tags$small("(Scroll down for instructions)")  )),
        solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
        p("Differential expression analysis is conducted for each hypothesis. There are two options to set (the adjusted-pvalue cut-off and the |logFC| cut-off).
          The results will appear in blocks (one per hypothesis) with 3 outputs:"),
        p("- the distribution of pvalue's: which has to be validated", a("(some help to identify the good shapes)", href="Pvalue_distrib.pdf"),""),
        p("- the MA plot (DE genes in red will varie with the p-value cutoff)"),
        p("- the table of statistics per gene/protein/metabolite (Number of stats displayed will varie with the p-value cutoff)")
    )
  })
  
  # diff param
  output$DiffParamUI <- renderUI({
    
    #we must run process before
    validate(
      need(rea.values[[dataset]]$process != FALSE, "Please run data processing")
    )
    
    validate(
      need(!is.null(rea.values$Contrasts.Sel), "Please run data processing ")
    )
    
    # #we must select list of contrast to test
    # validate(
    #   need(rea.values$analysis != FALSE, "Please select contrast")
    # )
    
    #design must be complete
    validate(
      need(rea.values[[dataset]]$compCheck != FALSE, session$userData$FlomicsMultiAssay@metadata$completeCheck[["error"]])
    )
    
    box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
        
        fluidRow(column(12,
                        ## list of contrasts to test
                        pickerInput(
                          inputId  = session$ns("contrastList"),
                          label    = "Selected contrast :",
                          choices  = session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                          multiple = TRUE, selected = session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName),
                        
                        # method for Diff analysis
                        selectInput(inputId  = session$ns("AnaDiffMethod"), label = "Method:",
                                    choices  = method,
                                    selected = method),
                        
                        numericInput(inputId = session$ns("Adj.pvalue.cutoff"), 
                                     label="Adjusted pvalue cutoff:",
                                     value=local.rea.values$Adj.pvalue.cutoff, min=0, max=1, 0.01),
                        numericInput(inputId = session$ns("abs.logFC.cutoff"),
                                     label="|logFC| cutoff:",
                                     value=local.rea.values$abs.logFC.cutoff, min=0, max=100, 0.01),
                        
                        # use of cluster. need setting step
                        materialSwitch(inputId = session$ns("clustermq"),
                                       label   =  popify(actionLink("infoCluster",paste0("Cluster: (?)")),"",
                                                         "If there is a huge number of contrasts, the calculation can be send to the cluster to be run in parrallel",
                                                         options=list(container="body"))
                                       , value = FALSE, status = "success"),
                        
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
    lapply(1:length(rea.values$Contrasts.Sel$contrast), function(i) {
      
      vect     <- unlist(rea.values$Contrasts.Sel[i,])
      
      # update/adapt PCA axis
      callModule(UpdateRadioButtons, paste0(vect["contrastName"],"-diff"))
      # select factors for color PCA plot
      #callModule(RadioButtonsCondition, paste0(vect["contrastName"],"-diff"), typeFact = c("Bio"))
    })
  })
  
  ##================================ RUN =======================================##
  
  ### run diff
  #######################
  observeEvent(input$runAnaDiff, {
    
    # check list of genes
    if(length(input$contrastList) == 0){
      
      showModal(modalDialog( title = "Error message", "Please select at least 1 hypothesis"))
    }
    validate({
      need(length(input$contrastList) != 0, message="Please select at least 1 hypothesis")
    })
    
    rea.values[[dataset]]$diffAnal   <- FALSE
    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    rea.values[[dataset]]$DiffValidContrast <- NULL
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]] <- NULL
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    if (dataset %in% rea.values$datasetDiff){
      rea.values$datasetDiff <- rea.values$datasetDiff[-which(rea.values$datasetDiff == dataset)]
    }
    
    
    dataset.SE@metadata$CoExpAnal   <- list()
    dataset.SE@metadata$DiffExpEnrichAnal  <- list()
    dataset.SE@metadata$CoExpEnrichAnal  <- list()
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    if(is.null(dataset.SE@metadata$DiffExpAnal)){
      
      dataset.SE@metadata$DiffExpAnal <- list()
      
      print(paste("# 4- Diff Analysis...", dataset))
      print(paste("#    => Filter Diff Analysis..."))
      
      # run diff analysis with selected method
      dataset.SE <- RunDiffAnalysis(object             = dataset.SE,
                                    design             = session$userData$FlomicsMultiAssay@metadata$design,
                                    contrastList       = session$userData$FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                                    Adj.pvalue.method  = "BH",
                                    DiffAnalysisMethod = input$AnaDiffMethod,
                                    clustermq          = input$clustermq,
                                    Adj.pvalue.cutoff  = input$Adj.pvalue.cutoff, 
                                    logFC.cutoff       = input$abs.logFC.cutoff)
    }
    else{
      
      print(paste("#    => Filter Diff Analysis..."))
      ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
      dataset.SE <- FilterDiffAnalysis(object            = dataset.SE,
                                       Adj.pvalue.cutoff = input$Adj.pvalue.cutoff,
                                       logFC.cutoff      = input$abs.logFC.cutoff)
    }
    
    # error management
    if(isFALSE(dataset.SE@metadata$DiffExpAnal[["results"]])){
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$DiffExpAnal[["ErrorStats"]])){
                               renderDataTable(dataset.SE@metadata$DiffExpAnal[["ErrorStats"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$DiffExpAnal[["Error"]])
                             }
      ))
    }
    
    if(is.null(dataset.SE@metadata$DiffExpAnal[["RawDEFres"]])){
      
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$DiffExpAnal[["ErrorTab"]])){
                               renderDataTable(dataset.SE@metadata$DiffExpAnal[["ErrorTab"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$DiffExpAnal[["error"]])
                             }
      ))
    }
    
    session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]] <- dataset.SE
    local.rea.values$dataset.SE <- dataset.SE
    
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
    
    
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]] 
    dataset.SE@metadata$DiffExpEnrichAnal <- NULL
    dataset.SE@metadata$CoExpEnrichAnal   <- NULL
    
    # filter DEG according pvalue adj cut-off
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal <-
      dataset.SE@metadata$DiffExpAnal
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]] <-
      rea.values[[dataset]]$DiffValidContrast
    
    rea.values[[dataset]]$diffValid <- TRUE
    rea.values$datasetDiff <- unique(c(rea.values$datasetDiff , dataset))
    
  }, ignoreInit = TRUE)
  
  ##============================== DISPLAY =====================================##
  
  # display results per contrast
  output$ContrastsResults <- renderUI({
    
    if (rea.values[[dataset]]$diffAnal == FALSE ||
        is.null(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["TopDEF"]])) return()
    
    # ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
    # local.rea.values$dataset.SE <- FilterDiffAnalysis(object = local.rea.values$dataset.SE,
    #                                                   Adj.pvalue.cutoff = input$Adj.pvalue.cutoff,
    #                                                   logFC.cutoff = input$abs.logFC.cutoff)
    # #Contrasts.Sel <- local.rea.values$dataset.SE@metadata$DiffExpAnal$contrasts
    
    # list of bio factors
    factors.bio   <- names(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("Bio")])
    factors.batch <- names(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("batch")])
    
    
    list(
      lapply(1:length(rea.values$Contrasts.Sel$contrast), function(i) {
        
        dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset, ".filtred")]]
        vect     <- unlist(rea.values$Contrasts.Sel[i,])
        res      <- dataset.SE@metadata$DiffExpAnal[["RawDEFres"]][[vect["contrastName"]]]
        stats    <- dataset.SE@metadata$DiffExpAnal[["stats"]][[vect["contrastName"]]]
        
        diff.plots <- DiffAnal.plot(dataset.SE, hypothesis=vect["contrastName"])
        
        if(dim(dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]])[1] == 0){
          
          fluidRow(
            column(10,
                   box(width=14, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "danger",
                       title = tags$h5(paste0(vect["tag"], " : ", vect["contrastName"],"  [#DE: ", stats$gDE," ]")),
                       
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
                       title = tags$h5(paste0(vect["tag"], ": ", vect["contrastName"],"  [#DE: ", stats$gDE," (up: ", stats$pgDEup,"%, ", "down: ", stats$pgDEdown,"%)]")),
                       
                       tabsetPanel(
                         
                         ### pvalue plot ###
                         tabPanel("Pvalue's distribution", renderPlot({ diff.plots$Pvalue.hist })),
                         
                         ### MAplot
                         tabPanel("MA plot", renderPlot({ suppressMessages(diff.plots$MA.plot) })),
                         
                         ### MAplot
                         tabPanel("Volcano plot", renderPlot({ suppressMessages(diff.plots$Volcano.plot) }, height = 600)),
                         
                         ### DEF result table ###
                         tabPanel("Table",
                                  ### DEF result table ###
                                  DT::renderDataTable({
                                    resTable <- dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]]
                                    resTable %>% DT::datatable(extensions = 'Buttons',
                                                               options = list(dom = 'lfrtipB',
                                                                              rownames = FALSE,
                                                                              pageLength = 10,
                                                                              buttons = c('csv', 'excel'),
                                                                              lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All")))) %>%
                                      formatStyle('logFC',
                                                  backgroundColor = styleInterval(c(0, 0.01), c('royalblue', 'white', 'red2')),
                                                  fontWeight = 'bold') %>% formatSignif(columns = 1:dim(resTable)[2], digits = 3)
                                  }, server = FALSE)),
                         
                         ### Heatmap ###
                         tabPanel("Heatmap",
                                  renderPlot({
                                    heatmap.plot(object=dataset.SE, hypothesis=vect["contrastName"], condition=input[[paste0(vect["contrastName"],"-heat.condColorSelect")]])
                                  }),
                                  renderText("Clustering method=ward.D2, center=TRUE, scale=FALSE"),
                                  ## select cluster to plot
                                  radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-heat.condColorSelect")),
                                               label = 'Levels :',
                                               choices = c("none", names(session$userData$FlomicsMultiAssay@colData)),
                                               selected = "none", inline = TRUE)
                         ),
                         
                         ### PCA ###
                         tabPanel("PCA on DE",
                                  
                                  fluidRow(
                                    column(width = 12,
                                           renderPlot({
                                             
                                             # 230317 :
                                             newDataset.SE <- dataset.SE[, which(colnames(dataset.SE) %in% dataset.SE@metadata$Groups$samples)]
                                             newDataset.SE <-  RunPCA(newDataset.SE[row.names(newDataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]])]) 
                                             
                                             PC1.value <- as.numeric(input[[paste0(vect["contrastName"],"-diff-Firstaxis")]][1])
                                             PC2.value <- as.numeric(input[[paste0(vect["contrastName"],"-diff-Secondaxis")]][1])
                                             condGroup <- input[[paste0(vect["contrastName"],"-pca.DE.condColorSelect")]][1]
                                             
                                              RFLOMICS::plotPCA(newDataset.SE,  PCA = "norm", PCs = c(PC1.value, PC2.value), condition = condGroup)
                                           })
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 6, 
                                           #RadioButtonsConditionUI(session$ns(paste0(vect["contrastName"],"-diff")))
                                           radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-pca.DE.condColorSelect")),
                                                        label = 'Levels:',
                                                        choices = c("groups",factors.bio),
                                                        selected = "groups")),
                                    column(width = 6, UpdateRadioButtonsUI(session$ns(paste0(vect["contrastName"],"-diff"))))
                                    
                                  )
                         ),
                         
                         ### boxplot DE ###
                         tabPanel("boxplot DE",
                                  fluidRow(
                                    column(width = 3,
                                           # selectInput(
                                           #   inputId = session$ns(paste0(vect["contrastName"], "-DE")), label = "Select DE:",
                                           #   choices = rownames(dplyr::arrange(dataset.SE@metadata$DiffExpAnal$TopDEF[[vect["contrastName"]]], Adj.pvalue)),
                                           #   multiple = FALSE, selectize = FALSE,
                                           #   size = 5 ),
                                           
                                           selectizeInput(
                                             inputId = session$ns(paste0(vect["contrastName"], "-DE")), label = "Select DE:",
                                             choices = rownames(dplyr::arrange(dataset.SE@metadata$DiffExpAnal$TopDEF[[vect["contrastName"]]], Adj.pvalue)),
                                             multiple = FALSE),
                                           
                                           #RadioButtonsConditionUI(session$ns(paste0(vect["contrastName"],"-DEcondition")))
                                           radioButtons(inputId = session$ns(paste0(vect["contrastName"],"-DEcondition")),
                                                        label = 'Levels:',
                                                        choices = c("groups",factors.bio),
                                                        selected = factors.bio[1])
                                           ),
                                    column(width = 9,
                                           renderPlot({
                                             boxplot.DE.plot(object=dataset.SE, DE=input[[paste0(vect["contrastName"], "-DE")]], 
                                                             condition=input[[paste0(vect["contrastName"], "-DEcondition")]]) })
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
    
    if (rea.values[[dataset]]$diffAnal == FALSE ||
        is.null(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["mergeDEF"]])) return()
    
    
    #local.rea.values$dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    DEF_mat <- as.data.frame(local.rea.values$dataset.SE@metadata$DiffExpAnal[["mergeDEF"]])
    
    index <- sapply(names(DEF_mat)[-1], function(x){(input[[paste0("checkbox_",x)]])}) %>% unlist()
    #index <- index[!sapply(index,is.null)]
    
    H_selected <- names(DEF_mat)[-1][index]
    
    rea.values[[dataset]]$DiffValidContrast <- dplyr::filter(rea.values$Contrasts.Sel, tag %in% H_selected)
    
    if (length(H_selected) > 1 && dim(DEF_mat)[1] != 0){
      
      box(width=12,  status = "warning",
          
          renderPlot({ UpSetR::upset(DEF_mat, sets = H_selected, order.by = "freq") })
      )
    }
  })
  
  return(input)
}



