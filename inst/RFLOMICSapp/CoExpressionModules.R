CoSeqAnalysisUI <- function(id){


  #name space for id
  ns <- NS(id)

  tagList(

    fluidRow(
      box(title = span(tagList(icon("cogs"), "   CoSeq ",a("(?)", href="https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html")  )), 
          solidHeader = TRUE, status = "warning", width = 12,
          
          "instructions... recommandation..."
          
      )),
    ### parametres for Co-Exp
    fluidRow(
      column(4, uiOutput(ns("CoExpParamUI"))),
      column(8, uiOutput(ns("CoExpResultUI"))))
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
               model <- "Normal"
               Trans <- "arcsin"
               normF <- "TMM"
               Gaussian <- "Gaussian_pk_Lk_Ck"},
           
           "metabolomics" = {
               model <- c("Normal","kmeans")
               Trans <- "none"
               normF <- "none"
               Gaussian <- c("Gaussian_pk_Lk_Ck", "Gaussian_pk_Lk_Bk", "none")},
           
           "proteomics" = {
               model <- c("Normal","kmeans")
               Trans <- "none"
               normF <- "none"
               Gaussian <- c("Gaussian_pk_Lk_Ck", "Gaussian_pk_Lk_Bk", "none")})
    
    names(model) <- paste0("model = ", model)
    names(Trans) <- paste0("transformation = ", Trans)
    names(normF) <- paste0("normFactors = ", normF)
    names(Gaussian) <- paste0("GaussianModel = ", Gaussian)
           
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
                               label    = "Gaussian Model (only for normal model) :", 
                               choices  = Gaussian, 
                               selected = Gaussian[1]))),
      
          fluidRow(
            column(8, 
                   sliderInput(session$ns("K.values"), label = "Number of clusters :", min=2, max=30, value=c(2,10), step=1)),
            
            column(4, 
                   numericInput(inputId = session$ns("iter"), label="Replicat :", value=10, 2, max=20, 1))),
      
   
          fluidRow(
            
            column(8,
                   materialSwitch(inputId = session$ns("clustermqCoseq"), label = "Cluster", value = FALSE, status = "success")),
                   # radioGroupButtons(inputId = session$ns("clustermqCoseq"), direction = "horizontal",
                   #                   label = " RUN :",
                   #                   choices = c("Local" = FALSE,  "genotoul" = TRUE),
                   #                   justified = FALSE, selected = "Local")
                   
                   
            column(4, actionButton(session$ns("runCoSeq"),"Run")))
    ))
  })
  
  
  # update K value (min max)
  observeEvent( input$K.values ,{
    
      min <- input$K.values[1]
      max <- input$K.values[2]
      if ((max - min) > 10) {max <- min + 10}
      
      
      # Control the value, min, max, and step.
      # Step size is 2 when input value is even; 1 when value is odd.
      updateSliderInput(session, "K.values", value = c(min, max),
                        min=2, max=30, step = 1)
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
                           as.character(dataset.SE@metadata$CoExpAnal[["error"]])))
    }

    #---- progress bar ----#
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#

    
    output$CoExpResultUI <- renderUI({
      
      # print coseq plots
      plot.coseq.res <- dataset.SE@metadata$CoExpAnal[["plots"]]
      
      box(title = "run clustering", status = "warning", solidHeader = TRUE, width = 14,
          
          tabBox( id = "runClustering", width = 12,
                  
                  tabPanel("logLike",  renderPlot({ plot.coseq.res$logLike })),
                  tabPanel("ICL",      renderPlot({ plot.coseq.res$ICL })),
                  tabPanel("profiles", renderPlot({ plot.coseq.res$profiles })),
                  tabPanel("boxplots", renderPlot({ plot.coseq.res$boxplots })),
                  tabPanel("boxplots_bis", 
                           renderUI({
                             
                             nb_cluster <- dataset.SE@metadata$CoExpAnal[["cluster.nb"]]
                             coseq.res  <- dataset.SE@metadata$CoExpAnal[["coseqResults"]]
                             
                             fluidRow(
                                 ## plot selected cluster(s)
                                 renderPlot({ coseq.y_profile.one.plot(coseq.res, input$selectCluster, dataset.SE@metadata$Groups) }),
                                 
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




