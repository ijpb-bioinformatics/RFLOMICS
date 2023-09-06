CoSeqAnalysisUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      box(title = span(tagList(icon('chart-line'), "   ",a("CoSeq", href="https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html"), tags$small("(Scroll down for instructions)")  )),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
          div(
            h4(tags$span("Parameters set up:", style = "color:orange")),
            p("You have first to choose between the ",tags$b("union")," or ",tags$b("intersection")," of your contrasts lists according to your biological question."),
            p("All the default parameters have been expertised according to each omic."),
            p("It is then recommanded to do a ",tags$b("first run")," with a large number of K with few iterations.
            If there is a K (Kbest different from Kmin and Kmax) for which the ICL is minimum (check the first graph obtained),
            then a second run has to be done with a larger experiment: the window of K can be centered around the Kbest and the number
            of technical replicates has to be increased to at least 20 iterations. For this larger experiment, it is recommanded to send
              analysis to remote ressources (see cluster option)"),
            h4(tags$span("Successful analysis:", style = "color:orange")),
            p("For a given K, the result will be considered as successful when at least the half of the iteration
              have run. Co-expression analysis will be considered as successful if there is at least a result for more than the
              half of K. In case of unsuccessful results, a detailed table of errors will appear."),
            h4(tags$span("Cluster option", style = "color:orange")),
            p("If you have a cluster account, you can configure a remote access to it
              (", a("see config_file", href="install_clustermq.txt"),")","and speed up results obtention."),
          ))),
    ### parametres for Co-Exp
    fluidRow(
      column(3, uiOutput(ns("CoExpParamUI"))),
      column(9, uiOutput(ns("CoExpResultUI")),
             uiOutput(ns("ERROR.UI"))))
  )
}

# tags$a(href="www.rstudio.com", "Click here!")

CoSeqAnalysis <- function(input, output, session, dataset, rea.values){
  
  local.rea.values <- reactiveValues(dataset.SE = NULL)
  
  # co-expression parameters
  output$CoExpParamUI <- renderUI({
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis")
    )
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    ##-> retrieve DEG lists and DEG valid lists
    ListNames.diff        <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]]$tag
    names(ListNames.diff) <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName
    
    ##-> option
    switch(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$omicType,
           
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
             Gaussian <- c("Gaussian_pk_Lk_Bk","Gaussian_pk_Lk_Ck",  "none")
             Scale <- TRUE
           },
           
           "proteomics" = {
             warning <- "(warning)"
             name <- "protein"
             model <- c("Normal","kmeans")
             Trans <- "none"
             normF <- "none"
             Gaussian <- c("Gaussian_pk_Lk_Bk", "Gaussian_pk_Lk_Ck","none")
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
      box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), status = "warning", width = 14,
          
          # Select lists of DGE to co-expression analysis
          fluidRow(
            column(12,
                   
                   pickerInput(
                     inputId  = session$ns("select"),
                     label    = "Validated DEG lists:",
                     choices  = ListNames.diff,
                     options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                     multiple = TRUE,
                     selected = ListNames.diff))),
          
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
                               label    = popify(actionLink("infoScale",paste0("Scale by ", name , ": (help)")),"","By default for proteomics or metabolomics data,coseq is done onto Z-scores (data scaled by proteins or metabolites) to group them according to their expression profile rather than abundance",options=list(container="body")),
                               choices  = Scale ,
                               selected = Scale[1])
            )
          ),
          
          fluidRow(
            column(12,
                   selectInput(session$ns("model"),
                               label    = "Default parameters:",
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
                               
                               label    = popify(actionLink("infoGaussianModel",paste0("Gaussian Model:", warning,"")),"","For proteomics or metabolomics data, coseq analysis may fail with default GaussianModel parameter. In this case, an error message will indicate to switch the other option: Gaussian_pk_Lk_Bk",options=list(container="body")),
                               choices  = Gaussian,
                               selected = Gaussian[1]))
          ),
          fluidRow(
            column(8,
                   sliderInput(session$ns("K.values"), label = "Number of clusters:", min=2, max=30, value=c(2,7), step=1)),
            column(4,
                   numericInput(inputId = session$ns("iter"), label="Iteration:", value=5, min = 5, max=20, step = 5))
          ),
          hr(),
          fluidRow(
            
            column(8,
                   materialSwitch(inputId = session$ns("clustermqCoseq"), label = "use remote cluster", value = FALSE, status = "success")
            ),
            column(4, actionButton(session$ns("runCoSeq"),"Run")))
      ))
  })
  
  #get list of DGE to process
  DEG_list <- reactive({
    getDEGlist_for_coseqAnalysis( matrix   = session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["mergeDEF"]],
                                  colnames = input$select, mergeType = input$unionInter)})
  
  # display nbr of selected genes
  output$mergeValue <- renderText({ 
    
    if(rea.values[[dataset]]$diffValid == FALSE) return()
    
    print(paste(length(DEG_list()), omics.dic[[session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$omicType]]$variableName, sep =" "))
  })
  
  # update K value (min max)
  observeEvent( c(input$K.values, input$clustermqCoseq, input$iter) ,{
    
    min <- input$K.values[1]
    max <- input$K.values[2]
    
    if(input$clustermqCoseq == FALSE && input$iter > 5){
      if ((max - min) > 5) {max <- min + 5}
    }
    
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session, "K.values", value = c(min, max),
                      min=2, max=30, step = 1)
    
  })
  
  ##================================ RUN =======================================##
  # run coexpression analysis
  # coseq
  observeEvent(input$runCoSeq, {
    
    print(paste0("# 9- CoExpression analysis... ", dataset ))
    
    rea.values[[dataset]]$coExpAnal  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    
    local.rea.values$dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    # check if no selected DGE list
    if(length(input$select) == 0){
      
      showModal(modalDialog( title = "Error message", "Please select at least 1 DEG list."))
    }
    validate({
      need(length(input$select) != 0, message="Please select at least 1 DEG list")
    })
    
    local.rea.values$dataset.SE@metadata$CoExpAnal   <- list()
    local.rea.values$dataset.SE@metadata$CoExpEnrichAnal  <- list()
    #FlomicsMultiAssay <- resetFlomicsMultiAssay(object=FlomicsMultiAssay, results=c("CoExpAnal"))
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run coseq: ", value = 0)
    on.exit(progress$close())
    progress$inc(1/2, detail = "In progress")
    #----------------------#
    
    
    # run coseq
    local.rea.values$dataset.SE <- runCoExpression(object = local.rea.values$dataset.SE, 
                                                   merge = input$unionInter, nameList = input$select,
                                                   K = input$K.values[1]:input$K.values[2], replicates = input$iter,
                                                   model  = input$model, transformation = input$transfo, normFactors = input$norm,
                                                   GaussianModel = input$GaussianModel, clustermq = input$clustermqCoseq, 
                                                   cmd = TRUE, silent = TRUE)
    
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]] <- local.rea.values$dataset.SE
    
    # If an error occured
    if(isFALSE(local.rea.values$dataset.SE@metadata$CoExpAnal[["results"]])){
      
      showModal(modalDialog( title = "Error message", paste0("No results! ", as.character(local.rea.values$dataset.SE@metadata$CoExpAnal[["error"]]))))
    }
    
    validate(
      need(!isFALSE(local.rea.values$dataset.SE@metadata$CoExpAnal[["results"]]), 
           paste0("No results!", as.character(local.rea.values$dataset.SE@metadata$CoExpAnal[["error"]])))
    )
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]] <- local.rea.values$dataset.SE
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#
    
    rea.values[[dataset]]$coExpAnal  <- TRUE
    
  }, ignoreInit = TRUE)
  
  ##============================== DISPLAY =====================================##
  
  output$CoExpResultUI <- renderUI({
    
    if (rea.values[[dataset]]$coExpAnal == FALSE) return()
    
    # Extract results
    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    factors.bio <- names(session$userData$FlomicsMultiAssay@metadata$design@Factors.Type[session$userData$FlomicsMultiAssay@metadata$design@Factors.Type %in% c("Bio")])
    
    plot.coseq.res <- dataset.SE@metadata$CoExpAnal[["plots"]]
    nb_cluster     <- dataset.SE@metadata$CoExpAnal[["cluster.nb"]]
    coseq.res      <- dataset.SE@metadata$CoExpAnal[["coseqResults"]]
    cluster.comp   <- dataset.SE@metadata$CoExpAnal[["clusters"]]
    topDEF         <- dataset.SE@metadata$DiffExpAnal$mergeDEF
    
    # For each id, get its cluster
    tab.clusters <- as.data.frame(ifelse(coseq.res@allResults[[names(nb_cluster)]] > 0.5, 1,0))
    tab.clusters <-  rownames_to_column(tab.clusters,var="DEF")
    Cluster.tab <- pivot_longer(data=tab.clusters,cols=2:(dim(tab.clusters)[2]),names_to="C",values_to="does.belong")
    
    # For each id, get its FC for all Contrasts
    tmp <- lapply(1:length(dataset.SE@metadata$DiffExpAnal$TopDEF),function(x){
      tmp <- dataset.SE@metadata$DiffExpAnal$TopDEF[x]
      df  <- data.frame("id"=row.names(tmp[[1]]), "logFC"=tmp[[1]]$logFC)
      names(df)=c("id",paste0("logFC.",filter(dataset.SE@metadata$DiffExpAnal$Validcontrasts,
                                              contrastName== names(tmp))$tag))
      return(df)
    })
    FC.tmp <- tmp %>% purrr::reduce(full_join, by='id')
    FC.tmp.mat <- as.matrix(FC.tmp[,-1])
    dimnames(FC.tmp.mat)[[1]] <- FC.tmp$id
    
    box(title = paste0("Number of clusters: ", nb_cluster), status = "warning", solidHeader = TRUE, width = 14,
        
        tabBox( id = "runClustering", width = 12,
                
                tabPanel("ICL",      
                         renderText("Integrated Completed Likelihood (ICL) plotted versus number of clusters"),
                         hr(),
                         renderPlot({ plot.coseq.res$ICL })),
                tabPanel("probapost_barplots", 
                         renderText("number of observations with a maximum conditional probability greater than threshold (0.8) per cluster"),
                         hr(),
                         renderPlot({ plot.coseq.res$probapost_barplots })),
                tabPanel("profiles",
                         renderPlot({ plot.coseq.res$boxplots + theme(axis.text.x = element_text(angle = 90, hjust = 1)) })),
                tabPanel("profiles per cluster",
                         fluidRow(
                           column(3, 
                                  ## select cluster to plot
                                  box(width = 14, background = "light-blue", title = NULL,
                                      radioButtons(inputId = session$ns("selectCluster"), label = "Select cluster:",
                                                   choices  = 1:nb_cluster, selected = 1),
                                      radioButtons(inputId = session$ns("profile.condition"), label = "Condition:",
                                                   choices  = c("groups", factors.bio), selected = "groups"),
                                      # selectInput(inputId = session$ns("observations"), label = "Observations(prob)",
                                      #             choices = cluster.comp[[input$selectCluster]], multiple = FALSE, selectize = FALSE, size = 5))
                                      uiOutput(session$ns("observationsUI"))
                                  )
                           ),
                           column(9, 
                                  renderPlot({ 
                                    coseq.profile.plot(dataset.SE, input$selectCluster, condition=input$profile.condition, observation=input$observations) }))
                         ),
                         hr(),
                         fluidRow(
                           renderText("blabla"),
                           renderPlot({ 
                             
                             tag <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]]$tag
                             
                             
                             if(length(tag) > 1){
                               UpSetR::upset(topDEF, sets = tag, order.by = "freq") 
                             }
                           })
                         )
                )
        )
    )
  })
  
  output$observationsUI  <- renderUI({

    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    coseq.res  <- dataset.SE@metadata$CoExpAnal[["coseqResults"]]
    clustr_num <- paste0("Cluster_",input$selectCluster)
    assays.data <- dplyr::filter(as.data.frame(coseq.res@assays@data[[1]]), get(clustr_num) > 0.8)# %>% dplyr::select(clustr_num)

    choices <- rownames(assays.data)
    names(choices) <- paste0(choices, " (",assays.data[,clustr_num], ")")
    
    selectizeInput(
      inputId = session$ns("observations"), label = "Observations(prob)",
      choices = choices, multiple = FALSE)
  })
  
  ## summary
  output$ERROR.UI <- renderUI({
    
    dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    if(rea.values[[dataset]]$coExpAnal == FALSE || dim(local.rea.values$dataset.SE@metadata$CoExpAnal[["stats"]])[1] == 0) return()
    
    box(title = "Failed cases", width = 14, status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
        
        DT::renderDataTable(DT::datatable(as.data.frame(dataset.SE@metadata$CoExpAnal[["stats"]]), 
                                          options = list(dom = 'tip'), rownames = FALSE))
    )
  })
}




