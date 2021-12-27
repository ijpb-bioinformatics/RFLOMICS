

DiffExpAnalysisUI <- function(id){

  #name space for id
  ns <- NS(id)

  tagList(
    fluidRow(
      uiOutput(ns("instruction"))),

    ### parametres for Diff Analysis
    fluidRow(
      column(3,uiOutput(ns("DiffParamUI")),
               uiOutput(ns("FilterPvalueUI"))),
      column(9,uiOutput(ns("ContrastsResults")),
               tags$br(),
               uiOutput(ns("ResultsMerge"))))

  )
}



DiffExpAnalysis <- function(input, output, session, dataset, rea.values){

  local.rea.values <- reactiveValues(dataset.SE = NULL)

  # list of tools for diff analysis
  MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")

  method <- switch (FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType,
                    "RNAseq"       = MethodList[1],
                    "proteomics"   = MethodList[2],
                    "metabolomics" = MethodList[2])

  output$instruction <- renderUI({
    box(title = span(tagList(icon("cogs"), "  ",  a(names(method), href="https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf"), "    (Scroll down for instructions)"  )),
        solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
        p("Differential expression analysis is conducted for each hypothesis. There is just two options to set (the ajusted-pvalue cut-off and the |logFC| cut-off).
          The results will appear in blocks (one per hypothesis) with 3 outputs:"),
        p("- the distribution of pvalue's : which has to be validated", a("(some help to identify the good shapes)", href="Pvalue_distrib.pdf"),""),
        p("- the MA plot (DE genes in red will varie with the p-value cutoff)"),
        p("- the table of statistics per gene/protein/metabolite (Number of stats displayed will varie with the p-value cutoff)")
    )
  })

  output$DiffParamUI <- renderUI({

    #we must select list of contrast to test
    validate(
      need(rea.values$analysis != FALSE, "Please select contrast")
    )

    box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",

        fluidRow(column(12,
               ## list of contrasts to test
               pickerInput(
                   inputId  = session$ns("contrastList"),
                   label    = "Selected contrast :",
                   choices  = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                   multiple = TRUE, selected = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName ),

               # method for Diff analysis
               selectInput(inputId  = session$ns("AnaDiffMethod"), label = "Method :",
                           choices  = method,
                           selected = method),

               # use of cluster. need setting step
               materialSwitch(inputId = session$ns("clustermq"),
                              label   =  popify(actionLink("infoCluster",paste0("Cluster: (?)")),"",
                                                "If there is a huge number of contrasts, the calculation can be send to the cluster to be run in parrallel",
                                                options=list(container="body"))
                              , value = FALSE, status = "success"),

               actionButton(session$ns("runAnaDiff"),"Run"))
    ))

  })

  # Run the differential analysis for each contrast set
  # Filter
  #   -> return a dynamic user interface with a collapsible box for each contrast
  #         - Pvalue graph
  #         - MAplot
  #         - Table of the DE genes
  #   -> combine data : union or intersection
  observeEvent(input$runAnaDiff, {

    # check list of genes
    if(length(input$contrastList) == 0){

      showModal(modalDialog( title = "Error message", "Please select at least 1 hypothesis"))
    }
    validate({
      need(length(input$contrastList) != 0, message="Please select at least 1 hypothesis")
    })

    print(paste("# 9- Diff Analysis...", dataset))

    rea.values[[dataset]]$diffAnal   <- FALSE
    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE

    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal <<- list()
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$CoExpAnal   <<- list()
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpEnrichAnal  <<- list()
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$CoExpEnrichAnal  <<- list()

    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#



    # run diff analysis with selected method
    local.rea.values$dataset.SE <- RunDiffAnalysis(object =             FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]],
                                                   design =             FlomicsMultiAssay@metadata$design,
                                                   contrastList =       FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                                                   Adj.pvalue.method =  "BH",
                                                   DiffAnalysisMethod = input$AnaDiffMethod,
                                                   clustermq =          input$clustermq)

    # error management
    if(isFALSE(local.rea.values$dataset.SE@metadata$DiffExpAnal[["results"]])){
      showModal(modalDialog( title = "Error message",
                             if(! is.null(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorStats"]])){
                               renderDataTable(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorStats"]],rownames = FALSE)
                             }
                             else{
                               as.character(local.rea.values$dataset.SE@metadata$DiffExpAnal[["Error"]])
                             }
      ))
    }

    if(is.null(local.rea.values$dataset.SE@metadata$DiffExpAnal[["RawDEFres"]])){

      showModal(modalDialog( title = "Error message",
                             if(! is.null(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorTab"]])){
                               renderDataTable(local.rea.values$dataset.SE@metadata$DiffExpAnal[["ErrorTab"]],rownames = FALSE)
                             }
                             else{
                               as.character(local.rea.values$dataset.SE@metadata$DiffExpAnal[["error"]])
                             }
      ))
    }

    rea.values[[dataset]]$diffAnal <- TRUE

    output$FilterPvalueUI <- renderUI({

      if (rea.values[[dataset]]$diffAnal == FALSE) return()

      box(title = NULL,status = "warning", width = 14,
          numericInput(inputId = session$ns("Adj.pvalue.cutoff"),
                   label="Adjusted pvalue cutoff :",
                   value=0.05, 0, max=1, 0.01),
          numericInput(inputId = session$ns("abs.logFC.cutoff"),
                       label="|logFC| cutoff :",
                       value=0, 0, max=1000, 0.1),
          actionButton(session$ns("validContrast"),"Validate"))
    })

    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#

  }, ignoreInit = TRUE)

  # display results per contrast
  output$ContrastsResults <- renderUI({

    if (rea.values[[dataset]]$diffAnal == FALSE) return()

    ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
    local.rea.values$dataset.SE <- FilterDiffAnalysis(object = local.rea.values$dataset.SE,
                                                      Adj.pvalue.cutoff = input$Adj.pvalue.cutoff,
                                                      logFC.cutoff = input$abs.logFC.cutoff)

    #Contrasts.Sel <<- local.rea.values$dataset.SE@metadata$DiffExpAnal$contrasts

    list(
      lapply(1:length(rea.values$Contrasts.Sel$contrast), function(i) {

      vect     <- unlist(rea.values$Contrasts.Sel[i,])
      res      <- local.rea.values$dataset.SE@metadata$DiffExpAnal[["RawDEFres"]][[vect["contrastName"]]]
      stats    <- local.rea.values$dataset.SE@metadata$DiffExpAnal[["stats"]][[vect["contrastName"]]]

      diff.plots <- DiffAnal.plot(local.rea.values$dataset.SE, hypothesis=vect["contrastName"],
                                  Adj.pvalue.cutoff = input$Adj.pvalue.cutoff, logFC.cutoff = input$abs.logFC.cutoff)

      fluidRow(
        column(10,
               box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning",
                   title = tags$h5(paste0(vect["tag"], " : ", vect["contrastName"],"  [#DE: ", stats$gDE," (up: ", stats$pgDEup,"%, ", "down: ", stats$pgDEdown,"%)]")),

                   tabsetPanel(

                    ### pvalue plot ###
                    tabPanel("Pvalue's distribution", renderPlot({ diff.plots$Pvalue.hist })),

                    ### MAplot
                    tabPanel("MA plot", renderPlot({ diff.plots$MA.plot })),

                    ### DEF result table ###
                    tabPanel("Table",
                       ### DEF result table ###
                       DT::renderDataTable({
                         resTable <- local.rea.values$dataset.SE@metadata$DiffExpAnal[["DEF"]][[vect["contrastName"]]]
                         keep <- resTable$Adj.pvalue <= input$Adj.pvalue.cutoff & abs(resTable$logFC) >= input$abs.logFC.cutoff
                         round(resTable[keep,],5) %>% DT::datatable(extensions = 'Buttons',
                                        options = list(dom = 'Blfrtip',
                                                      rownames = FALSE,
                                                      pageLength = 10,
                                                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                      lengthMenu = list(c(10,25,50,-1),c(10,25,50,"All")))) %>%
                           formatStyle('logFC',
                                 backgroundColor = styleInterval(c(0, 0.01), c('green', 'white', 'red')),
                                 fontWeight = 'bold')
                       })

                     )),
               )
        ),
        column(2,
               checkboxInput(session$ns(paste0("checkbox_", vect[["tag"]])), "OK", value = TRUE))
      )
    })
    )

  })

  # merge results on upset plot
  output$ResultsMerge <- renderUI({

    if (rea.values[[dataset]]$diffAnal == FALSE) return()


    index <- sapply(rea.values$Contrasts.Sel$tag, function(x){(input[[paste0("checkbox_",x)]])}) %>% unlist()

    H_selected <- rea.values$Contrasts.Sel$tag[index]

    DEF_mat <- as.data.frame(local.rea.values$dataset.SE@metadata$DiffExpAnal[["mergeDEF"]])

    rea.values[[dataset]]$DiffValidContrast <- dplyr::filter(rea.values$Contrasts.Sel, tag %in% H_selected)

    # Validcontrasts <- dplyr::filter(FlomicsMultiAssay@metadata$design@Contrasts.Sel, tag %in% H_selected)

    if (length(H_selected) > 1){

       box(width=12,  status = "warning",

           renderPlot({ UpSetR::upset(DEF_mat, sets = H_selected) })
       )
    }


  })

  # validate contrasts
  observeEvent(input$validContrast, {

    print(paste("# 9bis- Filter Diff Analysis...", dataset))

    rea.values[[dataset]]$diffValid  <- FALSE
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE

    # filter DEG according pvalue adj cut-off
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- local.rea.values$dataset.SE

    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]] <<- rea.values[[dataset]]$DiffValidContrast

    rea.values[[dataset]]$diffValid <- TRUE

  })

  return(input)
}



