CoSeqAnalysisUI <- function(id){


  #name space for id
  ns <- NS(id)

  tagList(

    ### parametres for Co-Exp
    fluidRow(
      column(4,
        box(title = "data", solidHeader = TRUE, status = "warning", width = 14,
            fluidRow(
              column(12,
                uiOutput(ns("selectDEGtoCoExp")))
            ),

            # Select type of merge : union or intersection
            fluidRow(
              column(12,
                selectInput(inputId = ns("unionInter"), label = "union or intersection",
                            choices = c("union", "intersection"), selected = "union"),
                verbatimTextOutput(ns("mergeValue"))
              )
            )
          ),
        box(title = span(tagList(icon("cogs"), "   CoSeq")), solidHeader = TRUE, status = "warning", width = 14,

          column(6,
            selectInput(ns("model"), label = "model :",
                        choices = list("normal"="normal","kmeans"="kmeans"), selected = "normal"),
            selectInput(ns("transfo"), label = "Transformation :",
                        choices = list("arcsin"="arcsin","none"="none"), selected = "arcsin"),
            selectInput(ns("norm"), label = "normFactors :",
                        choices = list("TMM"="TMM","none"="none"), selected = "TMM")
          ),
          column(6,
                 numericInput(inputId = ns("minK"), label="min K :", value=2 , 2, max=25, 1),
                 numericInput(inputId = ns("maxK"), label="max K :", value=20, 5, max=30, 1),
                 numericInput(inputId = ns("iter"), label="iteration :", value=5, 5, max=30, 1)
          )
        ),
        box(title = "run", solidHeader = TRUE, status = "warning", width = 14,
          column(6,
            selectInput(inputId = ns("clustermq-coseq"), label="send job to cluster",
                        choices = list("no"=FALSE,"genotoul"=TRUE))
          ),
          column(6, actionButton(ns("runCoSeq"),"Run clustering"))
        )
      ),
      column(8,
        box(title = "run clustering", status = "warning", solidHeader = TRUE, width = 14,

          tabBox( id = "runClustering", width = 12,

            tabPanel("logLike",  plotOutput(ns("logLike"))),
            tabPanel("ICL",      plotOutput(ns("ICL"))),
            tabPanel("profiles", plotOutput(ns("profiles"))),
            tabPanel("boxplots", plotOutput(ns("boxplots"))),
            tabPanel("boxplots_bis", uiOutput(ns("selectClusters"))),
            tabPanel("probapost_boxplots",  plotOutput(ns("probapost_boxplots"))),
            tabPanel("probapost_barplots",  plotOutput(ns("probapost_barplots"))),
            tabPanel("probapost_histogram", plotOutput(ns("probapost_histogram")))
          )
        )
      )
    )
  )
}

CoSeqAnalysis <- function(input, output, session, dataset){


  # Select lists of DGE to co-expression analysis
  output$selectDEGtoCoExp <- renderUI({

    ListNames.diff        <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["contrasts"]]$tag
    names(ListNames.diff) <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["contrasts"]]$contrastName

    pickerInput(
      inputId = session$ns("select"),
      label = "Select DEG lists:",
      choices = ListNames.diff,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
      multiple = TRUE, selected = ListNames.diff
    )
  })


  # get list of DGE to process
  # from union or intersection
  DEG_list <- reactive({getDEGlist_for_coseqAnalysis(matrix    = FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["mergeDEF"]],
                                           colnames  = input$select,
                                           mergeType = input$unionInter)})

  output$mergeValue <- renderText({ print(paste(length(DEG_list()), "genes", sep =" ")) })

  # run coexpression analysis
  # coseq
  observeEvent(input$runCoSeq, {

    # check if no selected DGE list
    if(length(input$select) == 0){

      showModal(modalDialog( title = "Error message", "Please select at least 1 DEG list."))
    }
    validate({
      need(length(input$select) != 0, message="Please select at least 1 DEG list")
    })


    print(paste("# 10- Co-expression analysis... ", dataset))

    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run coseq", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = paste("Doing part ", 10,"%", sep=""))
    #----------------------#

    # run coseq
    dataset.SE <- runCoExpression(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]],
                                  tools="coseq", geneList=DEG_list(),
                                  K=input$minK:input$maxK, iter = input$iter,
                                  model  = input$model,
                                  transformation=input$transfo, normFactors=input$norm,
                                  nameList=input$select, merge=input$unionInter)

    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- dataset.SE

    showModal(modalDialog( title = "Error message", dataset.SE@metadata$CoExpAnal[["error"]]))

    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#

    # print coseq plots
    plot.coseq.res <- dataset.SE@metadata$CoExpAnal[["plots"]]

    output$logLike  <- renderPlot({ plot.coseq.res$logLike })
    output$ICL      <- renderPlot({ plot.coseq.res$ICL })
    output$profiles <- renderPlot({ plot.coseq.res$profiles })
    output$boxplots <- renderPlot({ plot.coseq.res$boxplots })
    output$probapost_boxplots  <- renderPlot({ plot.coseq.res$probapost_boxplots })
    output$probapost_barplots  <- renderPlot({ plot.coseq.res$probapost_barplots })
    output$probapost_histogram <- renderPlot({ plot.coseq.res$probapost_histogram })

    #---- progress bar ----#
    progress$inc(3/4, detail = paste("Doing part ", 75,"%", sep=""))
    #----------------------#

    output$selectClusters <- renderUI({

      nb_cluster <- dataset.SE@metadata$CoExpAnal[["cluster.nb"]]
      coseq.res  <- dataset.SE@metadata$CoExpAnal[["coseqResults"]]

      fluidRow(

        ## plot selected cluster(s)
        renderPlot({ coseq.y_profile.one.plot(coseq.res, input$selectCluster, dataset.SE@metadata$Groups) }),

        ## select cluster to plot
        checkboxGroupInput(inputId = session$ns("selectCluster"), label = "Select cluster(s) :",
                           choices  = 1:nb_cluster, selected = 1, inline = TRUE)
      )
    })

    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#

  })
}
