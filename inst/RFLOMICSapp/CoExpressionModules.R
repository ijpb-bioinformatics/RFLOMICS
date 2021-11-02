CoSeqAnalysisUI <- function(id){


  #name space for id
  ns <- NS(id)

  tagList(

    fluidRow(
      box(title = span(tagList(icon("cogs"), "   CoSeq ",a("(see the guide)", href="https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html")  )),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = FALSE,
          div(
            h4(tags$span("Parameters set up:", style = "color:orange")),
            p("You have first to choose between the ",tags$b("union")," or ",tags$b("intersection")," of your contrasts lists according to your biological question."),
            p("All the default parameters have been expertised according to each omic."),
            p("It is then recommanded to do a ",tags$b("first run")," with a large number of K with few replicates.
            If there is a K (Kbest different from Kmin and Kmax) for which the ICL is minimum (check the first graph obtained),
            then a second run has to be done with a larger experiment: the window of K can be centered around the Kbest and the number
            of technical replicates has to be increased to at least 20 replicates. For this larger experiment, it is recommanded to send
              analysis to remote ressources (see cluster option)"),
            h4(tags$span("Successful analysis:", style = "color:orange")),
            p("For a given K, the result will be considered as successful when at least the half of the replicates
              have run. Co-expression analysis will be considered as successful if there is at least a result for more than the
              half of K. In case of unsuccessful results, a detailed table of errors will appear."),
            h4(tags$span("Cluster option", style = "color:orange")),
            p("If you have a cluster account, you can configure a remote access to it
              (", a("see config_file", href="install_clustermq.txt"),")","and speed up results obtention. Running coseq in local
              is very time consuming for the moment as iterations (one K, one replicate) run one after
              the other (ie. several hours for 5000 genes, a range of 20 K and 10 replicates),
              Whereas iteration run in parallel onto the cluster (ie, 2 minutes for the same manip).
              In this case, K and the number of genes have more influence on time"),
      ))),
    ### parametres for Co-Exp
    fluidRow(
      column(3, uiOutput(ns("CoExpParamUI"))),
      column(9, uiOutput(ns("CoExpResultUI"))))
  )
}

# tags$a(href="www.rstudio.com", "Click here!")

CoSeqAnalysis <- function(input, output, session, dataset){

  # co-expression parameters
  output$CoExpParamUI <- renderUI({

    # get good param :
    dataset.SE <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]

    ##-> retrieve DEG lists and DEG valid lists
    ListValidNames.diff   <- dataset.SE@metadata$DiffExpAnal[["Validcontrasts"]]$tag
    ListNames.diff        <- dataset.SE@metadata$DiffExpAnal[["contrasts"]]$tag
    names(ListNames.diff) <- dataset.SE@metadata$DiffExpAnal[["contrasts"]]$contrastName

    ##-> option
    switch(dataset.SE@metadata$omicType,

           "RNAseq" = {
               warning <- ""
               name <- "gene"
               model <- "Normal"
               Trans <- "arcsin"
               normF <- "TMM"
               Gaussian <- "Gaussian_pk_Lk_Ck"
               Scale <- FALSE
               },

           "metabolomics" = {
               warning <- "(warning)"
               name <- "metabolite"
               model <- c("Normal","kmeans")
               Trans <- "none"
               normF <- "none"
               Gaussian <- c("Gaussian_pk_Lk_Ck", "Gaussian_pk_Lk_Bk", "none")
               Scale <- TRUE
               },

           "proteomics" = {
               warning <- "(warning)"
               name <- "protein"
               model <- c("Normal","kmeans")
               Trans <- "none"
               normF <- "none"
               Gaussian <- c("Gaussian_pk_Lk_Ck", "Gaussian_pk_Lk_Bk", "none")
               Scale <- TRUE
               })

    names(model) <- paste0("model = ", model)
    names(Trans) <- paste0("transformation = ", Trans)
    names(normF) <- paste0("normFactors = ", normF)
    names(Gaussian) <- paste0("GaussianModel = ", Gaussian)
    names(Scale) <- paste0("Scale = ", Scale)

    # set param in interface
    tagList(

      ## Input parameters
      box(title = "Input", solidHeader = TRUE, status = "warning", width = 14,

          # Select lists of DGE to co-expression analysis
          fluidRow(
            column(12,

                   pickerInput(
                     inputId  = session$ns("select"),
                     label    = "Validated DEG lists:",
                     choices  = ListNames.diff,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = ListValidNames.diff))),

          # Select type of merge : union or intersection
          fluidRow(
            column(4,

                   radioButtons(inputId = session$ns("unionInter"), label=NULL ,
                                choices = c("union","intersection"),
                                selected = "union", inline = FALSE, width = 2)),

            column(8,
                   verbatimTextOutput(session$ns("mergeValue")) )),

          fluidRow(

            column(12,
                   selectInput(session$ns("scale"),
                               label    = popify(actionLink("infoScale",paste0("Scale by ", name , " : (help)")),"","By default for proteomics or metabolomics data,coseq is done onto Z-scores (data scaled by proteins or metabolites) to group them according to their expression profile rather than abundance",options=list(container="body")),
                               choices  = Scale ,
                               selected = Scale[1])
                   )
          ),

          fluidRow(
            column(12,
                   selectInput(session$ns("model"),
                               label    = "Default parameters :",
                               choices  = model ,
                               selected = model[1]),

                   selectInput(session$ns("transfo"),
                               label    = NULL,
                               choices  = Trans,
                               selected = Trans[1]),

                   selectInput(session$ns("norm"),
                               label    = NULL,
                               choices  = normF,
                               selected = normF[1]))),

          fluidRow(
            column(12,
                   selectInput(session$ns("GaussianModel"),

                               label    = popify(actionLink("infoGaussianModel",paste0("Gaussian Model :", warning,"")),"","For proteomics or metabolomics data, coseq analysis may fail with default GaussianModel parameter. In this case, an error message will indicate to switch the other option: Gaussian_pk_Lk_Bk",options=list(container="body")),
                               choices  = Gaussian,
                               selected = Gaussian[1]))),

          fluidRow(
            column(8,
                   sliderInput(session$ns("K.values"), label = "Number of clusters :", min=2, max=30, value=c(2,10), step=1)),
            column(4,
                   numericInput(inputId = session$ns("iter"), label="Replicat :", value=10, 2, max=20, 1))),


          fluidRow(

            column(8,
                   materialSwitch(inputId = session$ns("clustermqCoseq"), label = "Cluster", value = FALSE, status = "success"),
                   ),
                   # radioGroupButtons(inputId = session$ns("clustermqCoseq"), direction = "horizontal",
                   #                   label = " RUN :",
                   #                   choices = c("Local" = FALSE,  "genotoul" = TRUE),
                   #                   justified = FALSE, selected = "Local")


            column(4, actionButton(session$ns("runCoSeq"),"Run")))
    ))
  })


  # update K value (min max)
  observeEvent( input$K.values ,{

    switch(as.character(input$clustermqCoseq),{
      `FALSE`={
      min <- input$K.values[1]
      max <- input$K.values[2]
      if ((max - min) > 10) {max <- min + 10}

      # Control the value, min, max, and step.
      # Step size is 2 when input value is even; 1 when value is odd.
      updateSliderInput(session, "K.values", value = c(min, max),
                        min=2, max=30, step = 1)
      }
      `TRUE`={
        updateSliderInput(session, "K.values", value = c(min, max),
                          min=2, max=30, step = 1)
      }
    }
  )
})



  # get list of DGE to process
  DEG_list <- reactive({
    getDEGlist_for_coseqAnalysis( matrix   = FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["mergeDEF"]],
                                  colnames = input$select, mergeType = input$unionInter)})

  # display nbr of selected genes
  output$mergeValue <- renderText({
    print(paste(length(DEG_list()), "genes", sep =" ")) })

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
    dataset.SE <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]

    dataset.SE <- runCoExpression(object=dataset.SE, geneList=DEG_list(), merge=input$unionInter, nameList=input$select,
                                  K=input$K.values[1]:input$K.values[2], replicates = input$iter,
                                  model  = input$model, transformation=input$transfo, normFactors=input$norm,
                                  GaussianModel =input$GaussianModel, clustermq = input$clustermqCoseq)

    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- dataset.SE

    # If an error occured
    if(isFALSE(dataset.SE@metadata$CoExpAnal[["results"]]))
    {
      showModal(modalDialog( title = "Error message",
                             if(! is.null(dataset.SE@metadata$CoExpAnal[["stats"]])){
                             renderDataTable(dataset.SE@metadata$CoExpAnal[["stats"]],rownames = FALSE)
                             }
                             else{
                               as.character(dataset.SE@metadata$CoExpAnal[["error"]])
                             }
                             ))
    }

    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#


    output$CoExpResultUI <- renderUI({

      # print coseq plots
      plot.coseq.res <- dataset.SE@metadata$CoExpAnal[["plots"]]

      box(title = "run clustering", status = "warning", solidHeader = TRUE, width = 14,

          tabBox( id = "runClustering", width = 12,

                  tabPanel("ICL",      renderPlot({ plot.coseq.res$ICL })),
                  tabPanel("logLike",  renderPlot({ plot.coseq.res$logLike })),
                  tabPanel("profiles", renderPlot({ plot.coseq.res$profiles })),
                  tabPanel("boxplots", renderPlot({ plot.coseq.res$boxplots })),
                  tabPanel("boxplots_bis",
                           renderUI({

                             nb_cluster <- dataset.SE@metadata$CoExpAnal[["cluster.nb"]]
                             coseq.res  <- dataset.SE@metadata$CoExpAnal[["coseqResults"]]

                             fluidRow(
                                 ## plot selected cluster(s)
                                 renderPlot({ coseq.y_profile.one.plot(coseq.res, input$selectCluster, dataset.SE@metadata$Groups) }),
                                 ## print datatable of metabolite
                                 DT::renderDataTable(DT::datatable(as.data.frame(coseq.res@allResults[[1]]) %>%
                                                   select(.,paste0("Cluster_",input$selectCluster)) %>%
                                                   dplyr::filter(.,get(paste0("Cluster_",input$selectCluster))==1)),
                                                   options = list(rownames = FALSE, pageLength = 10)),
                                 ## select cluster to plot
                                 checkboxGroupInput(inputId = session$ns("selectCluster"), label = "Select cluster(s) :",
                                                    choices  = 1:nb_cluster, selected = 1, inline = TRUE))})),

                  tabPanel("probapost_boxplots",  renderPlot({ plot.coseq.res$probapost_boxplots })),
                  tabPanel("probapost_barplots",  renderPlot({ plot.coseq.res$probapost_barplots })),
                  tabPanel("probapost_histogram", renderPlot({ plot.coseq.res$probapost_histogram }))
          )
      )
    })


    #---- progress bar ----#
    progress$inc(3/4, detail = paste("Doing part ", 75,"%", sep=""))
    #----------------------#


    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#

  }, ignoreInit = TRUE)
}




