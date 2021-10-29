

DiffExpAnalysisUI <- function(id){

  #name space for id
  ns <- NS(id)

  tagList(
    fluidRow(
      box(title = "Instructions:",
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = FALSE,
        p("Differential expression/abundance analysis is conducted for each hypothesis. There is just one option to set (the ajusted-pvalue cut-off, which is set to 5 % by default).
        The results will appear in blocks (one per hypothesis) with 3 outputs:"),
        p("- the distribution of pvalue's : which has to be validated", a("(some help to identify the good shapes)", href="Pvalue_distrib.pdf"),""),
        p("- the MA plot (DE genes in red will varie with the p-value cutoff)"),
        p("- the table of statistics per gene/protein/metabolite (Number of stats displayed will varie with the p-value cutoff)")
      )),

    ### parametres for Diff Analysis
    fluidRow(
      column(3,uiOutput(ns("DiffParamUI"))),
      column(9,uiOutput(ns("ContrastsResults")))),
    tags$br(),
    fluidRow(
      column(3,uiOutput(ns("FilterPvalueUI"))),
      column(9,uiOutput(ns("ResultsMerge"))))
  )
}



DiffExpAnalysis <- function(input, output, session, dataset){

  output$DiffParamUI <- renderUI({

    MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")

    method <- switch (FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType,
                      "RNAseq"       = MethodList[1],
                      "proteomics"   = MethodList[2],
                      "metabolomics" = MethodList[2])

    box(title = span(tagList(icon("cogs"), "   ", method)), width = 14, status = "warning",

        fluidRow(column(12,
               ## list of contrasts to test

               pickerInput(
                   inputId  = session$ns("contrastList"),
                   label    = "Selected contrast :",
                   choices  = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                   multiple = TRUE, selected = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName ),

               selectInput(inputId  = session$ns("AnaDiffMethod"), label = "Method :",
                           choices  = method,
                           selected = method),

               materialSwitch(inputId = session$ns("clustermq"), label =  popify(actionLink("infoCluster",paste0("Cluster: (?)")),"",
                                                                                 "If there is a huge number of contrasts, the calculation can be send to the cluster to be run in parrallel",options=list(container="body"))
                              , value = FALSE, status = "success"),
              #selectInput(inputId = session$ns("clustermq"),
              #label   ="send job to cluster",
              #choices = list("no"=FALSE,"genotoul"=TRUE)),

              actionButton(session$ns("runAnaDiff"),"Run"))
  ))

  })

  # list of selected contrast
  # output$contrastListUI <- renderUI({
  #   pickerInput(
  #     inputId = session$ns("contrastList"),
  #     label = "Selected contrast :",
  #     choices = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
  #     options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
  #     multiple = TRUE, selected = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName )
  # })


  # chose of methods based on data type
  output$AnaDiffMethodUI <- renderUI({

    MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "lmFit (limma)"="limmalmFit")

    method <- switch (FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType,
                      "RNAseq"       = MethodList[1],
                      "proteomics"   = MethodList[2],
                      "metabolomics" = MethodList[2])

    selectInput(session$ns("AnaDiffMethod"), label = "Method :",
                choices  = method,
                selected = method)
  })

  # Run the differential analysis for each contrast set
  # Filter
  #   -> return a dynamic user interface with a collapsible box for each contrast
  #         - Pvalue graph
  #         - MAplot
  #         - Table of the DE genes
  #   -> combine data : union or intersection
  observeEvent(input$runAnaDiff, {
    print(paste0("observeEvent",input$runAnaDiff))

    # check list of genes
    if(length(input$contrastList) == 0){

      showModal(modalDialog( title = "Error message", "Please select at least 1 hypothesis"))
    }
    validate({
      need(length(input$contrastList) != 0, message="Please select at least 1 hypothesis")
    })

    print(paste("# 9- Diff Analysis...", dataset))

    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = paste("Doing part ", 10,"%", sep=""))
    #----------------------#

    # run diff analysis with select method
    dataset.SE <- RunDiffAnalysis(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]], design = Design,
                                  contrastList = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                                  Adj.pvalue.method="BH",
                                  DiffAnalysisMethod=input$AnaDiffMethod,
                                  clustermq=input$clustermq)

    if(isFALSE(dataset.SE@metadata$DiffExpAnal[["results"]]))
    {
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$DiffExpAnal[["ErrorStats"]])){
                               renderDataTable(dataset.SE@metadata$DiffExpAnal[["ErrorStats"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$DiffExpAnal[["Error"]])
                             }
      ))
    }



    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- dataset.SE

    if(is.null(dataset.SE@metadata$DiffExpAnal[["RawDEFres"]]))
    {
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$DiffExpAnal[["ErrorTab"]])){
                               renderDataTable(dataset.SE@metadata$DiffExpAnal[["ErrorTab"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$CoExpAnal[["error"]])
                             }
      ))
    }

    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#

    output$FilterPvalueUI <- renderUI({
      box(title = NULL,status = "warning", width = 10,
          numericInput(inputId = session$ns("Adj.pvalue.cutoff"),
                   label="Adjusted pvalue cutoff :",
                   value=0.05, 0, max=1, 0.01))
    })

  }, ignoreInit = TRUE)


  observeEvent((input$Adj.pvalue.cutoff ), {

    print(paste("# 9bis- Filter Diff Analysis...", dataset))

    ### adj_pvalue filtering by calling the RundDiffAnalysis method without filtering
    dataset.SE <- FilterDiffAnalysis(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]],
                                     Adj.pvalue.cutoff = input$Adj.pvalue.cutoff)
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- dataset.SE

    Contrasts.Sel <- Design@Contrasts.Sel

    output$ContrastsResults <- renderUI({
      list(
      lapply(1:length(Contrasts.Sel$contrast), function(i) {

        vect     <- unlist(Contrasts.Sel[i,])
        res      <- dataset.SE@metadata$DiffExpAnal[["RawDEFres"]][[vect["contrastName"]]]
        stats    <- dataset.SE@metadata$DiffExpAnal[["stats"]][[vect["contrastName"]]]

        diff.plots <- DiffAnal.plot(dataset.SE, hypothesis=vect["contrastName"], Adj.pvalue.cutoff = input$Adj.pvalue.cutoff)

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
                         resTable <- dataset.SE@metadata$DiffExpAnal[["DEF"]][[vect["contrastName"]]]
                         DT::datatable(round(resTable[resTable$Adj.pvalue <= input$Adj.pvalue.cutoff,],5), options = list(rownames = FALSE, pageLength = 10))
                       })
                     )),
                 )
          ),
          column(2,
                 checkboxInput(session$ns(paste0("checkbox_", vect[["tag"]])), "OK", value = TRUE))
        )
      }),

        fluidRow(
          column(10),
          column(2,
                 actionButton(session$ns("validContrast"),"Validate"))
        ))

    })


    # merge results on upset plot
    output$ResultsMerge <- renderUI({

      index <- sapply(Contrasts.Sel$tag, function(x){(input[[paste0("checkbox_",x)]])}) %>% unlist()

      H_selected <- Contrasts.Sel$tag[index]

      DEF_mat <- as.data.frame(dataset.SE@metadata$DiffExpAnal[["mergeDEF"]])

      Validcontrasts <- dplyr::filter(Design@Contrasts.Sel, tag %in% H_selected)
      FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]] <<- Validcontrasts

      if (length(H_selected) > 1){
        fluidRow(
          column(10,
                 box(width=12,  status = "warning",

                     renderPlot({ UpSetR::upset(DEF_mat, sets = H_selected) })
                 )
          )
        )
      }


    })

  }, ignoreInit = TRUE)


  observeEvent(input$validContrast, {

  })

  return(input)
}



