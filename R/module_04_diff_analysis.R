
#' @importFrom UpSetR upset
#' @importFrom DT formatSignif datatable renderDataTable styleInterval formatStyle
#' @importFrom htmltools span tagList p div a h4 h5 hr tags br HTML
#' @importFrom shinyBS popify bsButton addPopover bsTooltip
#' @importFrom shinydashboard box tabBox updateTabItems menuItem menuItemOutput 
#' tabItem renderMenu tabItems sidebarMenu menuSubItem
#' @rawNamespace import(shiny, except = renderDataTable)
#' @importFrom shinyWidgets pickerInput materialSwitch
#' @importFrom colourpicker colourInput
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
          There are two options to set (the adjusted-pvalue cut-off and the |logFC| cut-off).
          The results will appear in blocks (one per contrast) with several outputs:"),
        p("- Distribution of pvalue's which has to be validate to trust results", a("(some help to identify the good shapes)", href="/www/Pvalue_distrib.pdf"),""),
        p("- MA plot"),
        p("- Volcano plot"),
        p("- Dataframe with statistical results"),
        p("- Heatmap"),
        p("- PCA")
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
      need(rea.values[[dataset]]$compCheck != FALSE, session$userData$FlomicsMultiAssay@metadata$completeCheck[["error"]])
    )
    
    local.rea.values$selectedContrasts <- 
      getExpressionContrast(session$userData$FlomicsMultiAssay[[dataset]], modelFormula = getModelFormula(session$userData$FlomicsMultiAssay)) %>% 
      purrr::reduce(rbind) %>% dplyr::filter(contrast %in% rea.values$Contrasts.Sel$contrast)
    
    
    ## getcontrast
    box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
        
        fluidRow(column(12,
                        
                        ## list of contrasts to test
                        pickerInput(
                          inputId  = session$ns("contrastList"),
                          label    = "Selected contrasts:",
                          choices  = local.rea.values$selectedContrasts$contrastName),
                        
                        # method for Diff analysis
                        selectInput(inputId  = session$ns("AnaDiffMethod"), label = "Method:",
                                    choices  = method,
                                    selected = method),
                        
                        numericInput(inputId = session$ns("p.adj.cutoff"), 
                                     label="Adjusted pvalue cutoff:",
                                     value=local.rea.values$p.adj.cutoff, min=0, max=1, 0.01),
                        numericInput(inputId = session$ns("abs.logFC.cutoff"),
                                     label="|logFC| cutoff:",
                                     value=local.rea.values$abs.logFC.cutoff, min=0, max=100, 0.01),
                        
                        # use of cluster. need setting step
                        materialSwitch(inputId = session$ns("clustermq"),
                                       label   =  popify(actionLink(inputId=session$ns("infoCluster"),label=paste0("Cluster: (?)")),
                                                         title= "",
                                                         content="If there is a huge number of contrasts, the calculation can be send to the cluster to be run in parrallel",
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
    rea.values[[dataset]]$DiffValidContrast <- NULL
    
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpEnrichAnal <- NULL
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpAnal         <- NULL
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpEnrichAnal   <- NULL
    
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal[["Validcontrasts"]] <- NULL
    
    if (dataset %in% rea.values$datasetDiff){
      rea.values$datasetDiff <- rea.values$datasetDiff[-which(rea.values$datasetDiff == dataset)]
    }
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    if(is.null(dataset.SE@metadata$DiffExpAnal)){
      
      dataset.SE@metadata$DiffExpAnal <- list()
      
      print(paste("# 4- Diff Analysis...", dataset))
      print(paste("#    => Filter Diff Analysis..."))
      
      # run diff analysis with selected method
      dataset.SE <- runDiffAnalysis(object             = dataset.SE,
                                    design             = session$userData$FlomicsMultiAssay@metadata$design,
                                    modelFormula       = session$userData$FlomicsMultiAssay@metadata$design$Model.formula,
                                    contrastList       = rea.values$Contrasts.Sel,
                                    p.adj.method  = "BH",
                                    method = input$AnaDiffMethod,
                                    clustermq          = input$clustermq,
                                    p.adj.cutoff  = input$p.adj.cutoff, 
                                    logFC.cutoff       = input$abs.logFC.cutoff,
                                    cmd                = TRUE)
    }
    else{
      
      print(paste("#    => Filter Diff Analysis..."))
      ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
      dataset.SE <- filterDiffAnalysis(object            = dataset.SE,
                                       p.adj.cutoff = input$p.adj.cutoff,
                                       logFC.cutoff      = input$abs.logFC.cutoff)
    }
    
    # error management
    if(isFALSE(dataset.SE@metadata$DiffExpAnal[["results"]])){
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$DiffExpAnal[["ErrorStats"]])){
                               DT::renderDataTable(dataset.SE@metadata$DiffExpAnal[["ErrorStats"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$DiffExpAnal[["Error"]])
                             }
      ))
    }
    
    if(is.null(dataset.SE@metadata$DiffExpAnal[["RawDEFres"]])){
      
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$DiffExpAnal[["ErrorTab"]])){
                               DT::renderDataTable(dataset.SE@metadata$DiffExpAnal[["ErrorTab"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$DiffExpAnal[["error"]])
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
    
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpEnrichAnal <- NULL
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpEnrichAnal   <- NULL
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$CoExpAnal         <- NULL
    
    session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal[["Validcontrasts"]] <- rea.values[[dataset]]$DiffValidContrast
    
    rea.values[[dataset]]$diffValid <- TRUE
    rea.values$datasetDiff <- unique(c(rea.values$datasetDiff , dataset))
    
  }, ignoreInit = TRUE)
  
  ##============================== DISPLAY =====================================##
  
  # display results per contrast
  output$ContrastsResults <- renderUI({
    
    if (rea.values[[dataset]]$diffAnal == FALSE ||
        is.null(session$userData$FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal[["TopDEF"]])) return()
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
    
    # list of bio factors
    factors.bio   <- bioFactors(dataset.SE)
    factors.batch <- batchFactors(dataset.SE)
    factors.meta  <- metaFactors(dataset.SE)
    
    list(
      lapply(1:nrow(getSelectedContrasts(dataset.SE)), function(i) {
        
        vect     <- unlist(dataset.SE@metadata$design$Contrasts.Sel[i,])
        res      <- dataset.SE@metadata$DiffExpAnal[["RawDEFres"]][[vect["contrastName"]]]
        stats    <- dataset.SE@metadata$DiffExpAnal[["stats"]][vect["contrastName"],]
        
        diff.plots <- plotDiffAnalysis(dataset.SE, contrastName=vect["contrastName"])
        
        if(dim(dataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]])[1] == 0){
          
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
                       title = tags$h5(paste0(vect["tag"], ": ", vect["contrastName"],"  [#DE: ", stats[["All"]]," ; Up: ", stats[["Up"]]," (",round(stats[["Up"]]/stats[["All"]],2)*100," %)"," ; ", "Down: ", stats[["Down"]]," (",round(stats[["Down"]]/stats[["All"]],2)*100,"%)")),
                       
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
                                  renderText("Clustering method = ward.D2, center = TRUE, scale = FALSE"),
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
                                  
                                  fluidRow(
                                    column(width = 12,
                                           renderPlot({
                                             
                                             # 230317 :
                                             newDataset.SE <- dataset.SE[, which(colnames(dataset.SE) %in% dataset.SE@metadata$Groups$samples)]
                                             newDataset.SE <-  runOmicsPCA(newDataset.SE[row.names(newDataset.SE@metadata$DiffExpAnal[["TopDEF"]][[vect["contrastName"]]])],ncomp = 5, raw = FALSE) 
                                             
                                             PC1.value <- as.numeric(input[[paste0(vect["contrastName"],"-diff-Firstaxis")]][1])
                                             PC2.value <- as.numeric(input[[paste0(vect["contrastName"],"-diff-Secondaxis")]][1])
                                             condGroup <- input[[paste0(vect["contrastName"],"-pca.DE.condColorSelect")]][1]
                                             
                                             RFLOMICS::plotPCA(newDataset.SE,  raw = "norm", axes = c(PC1.value, PC2.value), groupColor = condGroup)
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
                            
                                           selectizeInput(
                                             inputId = session$ns(paste0(vect["contrastName"], "-DE")), label = "Select DE:",
                                             choices = rownames(dplyr::arrange(dataset.SE@metadata$DiffExpAnal$TopDEF[[vect["contrastName"]]], Adj.pvalue)),
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
        is.null(dataset.SE@metadata$DiffExpAnal[["mergeDEF"]])) return()
    
    DEF_mat <- as.data.frame(dataset.SE@metadata$DiffExpAnal[["mergeDEF"]])
    
    index <- sapply(names(DEF_mat)[-1], function(x){(input[[paste0("checkbox_",x)]])}) %>% unlist()
    
    H_selected <- names(DEF_mat)[-1][index]
    
    rea.values[[dataset]]$DiffValidContrast <- dplyr::filter(dataset.SE@metadata$design$Contrasts.Sel, tag %in% H_selected)
    
    if (length(H_selected) > 1 && dim(DEF_mat)[1] != 0){
      
      box(width=12,  status = "warning",
          
          renderPlot({ UpSetR::upset(DEF_mat, sets = H_selected, order.by = "freq") })
      )
    }
  })
  
  return(input)
}

############## functions ###############

# ----- check run diff execution ------
check_run_diff_execution <- function(object.SE, param.list = NULL){
  
  # filtering setting
  if(is.null(object.SE@metadata[["DiffExpAnal"]]) || object.SE@metadata[["DiffExpAnal"]]$results != TRUE) return(TRUE)
  
  if(param.list$method             != getDiffSettings(object.SE)$method)            return(TRUE)
  if(param.list$p.adj.method  != getDiffSettings(object.SE)$p.adj.method) return(TRUE)
  if(param.list$p.adj.cutoff  != getDiffSettings(object.SE)$p.adj.cutoff) return(TRUE)
  if(param.list$abs.logFC.cutoff   != getDiffSettings(object.SE)$abs.logFC.cutoff)  return(TRUE)
  
  return(FALSE)
}

