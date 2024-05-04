CoSeqAnalysisUI <- function(id){
    
    #name space for id
    ns <- NS(id)
    
    tagList(
        
        fluidRow(
            box(title = span(
                tagList(icon('chart-line'), "   ",
                        a("CoSeq", href="https://www.bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html"), 
                        tags$small("(Scroll down for instructions)")  )),
                solidHeader = TRUE, status = "warning", width = 12, 
                collapsible = TRUE, collapsed = TRUE,
                div(
                    p("Analyses in this module are conducted using the 
                      coseq R-package. If you have more questions or
                      interest in this package, 
              please check the associated paper or the online vignette
              at https://bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html"),
              h4(tags$span("Parameters set up:", style = "color:orange")),
              p("You have first to choose between the ",tags$b("union"),
                " or ",tags$b("intersection")," of your contrasts lists 
                according to your biological question."),
              p("All the default parameters have been expertised according 
                to each omic."),
              p("It is then recommanded to do a ",tags$b("first run"),
                " with a large number of K with few iterations.
            If there is a K (Kbest different from Kmin and Kmax) for which 
            the ICL is minimum (check the first graph obtained),
            then a second run has to be done with a larger experiment: 
            the window of K can be centered around the Kbest and the number
            of technical replicates has to be increased to at least 20
            iterations. For this larger experiment, it is recommanded to send
              analysis to remote ressources (see cluster option)"),
            h4(tags$span("Successful analysis:", style = "color:orange")),
            p("For a given K, the result will be considered as successful
              when at least the half of the iteration
              have run. Co-expression analysis will be considered as 
              successful if there is at least a result for more than the
              half of K. In case of unsuccessful results, a detailed table 
              of errors will appear."),
            h4(tags$span("Cluster option", style = "color:orange")),
            p("If you have a cluster account, you can configure a remote 
              access to it
              (", a("see config_file", href="install_clustermq.txt"),")",
              "and speed up results obtention."),
                ))),
        ### parametres for Co-Exp
        fluidRow(
            column(3, uiOutput(ns("CoExpParamUI"))),
            column(9, uiOutput(ns("CoExpResultUI")),
                   uiOutput(ns("ERROR.UI"))))
    )
}

#' @importFrom tidyr pivot_longer 
#' @importFrom tibble rownames_to_column
#' @importFrom DT renderDataTable  datatable
#' @importFrom dplyr filter
#' @importFrom purrr reduce
#' @importFrom UpSetR upset
CoSeqAnalysis <- function(input, output, session, dataset, rea.values){
    
    local.rea.values <- reactiveValues(features.list = NULL)
    
    # co-expression parameters
    output$CoExpParamUI <- renderUI({
        
        validate(
            need(rea.values[[dataset]]$diffValid != FALSE, 
                 "Please run diff analysis")
        )
        
        if(rea.values[[dataset]]$diffValid == FALSE) return()
        
        MAE.data <- session$userData$FlomicsMultiAssay
        SE.data <- MAE.data[[dataset]]
        SE.filtered <- MAE.data[[dataset]]
        
        ##-> retrieve DEG lists and DEG valid lists
        ListNames.diff        <- getValidContrasts(SE.filtered)$contrastName
        names(ListNames.diff) <- getValidContrasts(SE.filtered)$contrastName
        
        ##-> option
        switch(getOmicsTypes(SE.filtered),
               
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
                   Gaussian <- c("Gaussian_pk_Lk_Bk",
                                 "Gaussian_pk_Lk_Ck",  "none")
                   Scale <- TRUE
               },
               
               "proteomics" = {
                   warning <- "(warning)"
                   name <- "protein"
                   model <- c("Normal","kmeans")
                   Trans <- "none"
                   normF <- "none"
                   Gaussian <- c("Gaussian_pk_Lk_Bk", 
                                 "Gaussian_pk_Lk_Ck","none")
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
            box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), 
                status = "warning", width = 14,
                
                # Select lists of DGE to co-expression analysis
                fluidRow(
                    column(12,
                           
                           pickerInput(
                               inputId  = session$ns("select"),
                               label    = .addBSpopify(label = 'Validated DEG lists:', 
                                                       content = paste0("Choose between the union or intersection",
                                                                        "of your contrasts lists according to your biological question.")),
                               choices  = ListNames.diff,
                               options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                               multiple = TRUE,
                               selected = ListNames.diff))),
                
                # Select type of merge : union or intersection
                fluidRow(
                    column(4,
                           
                           radioButtons(inputId = session$ns("unionInter"), 
                                        label=NULL ,
                                        choices = c("union","intersection"),
                                        selected = "union", inline = FALSE, width = 2)),
                    
                    column(8,
                           verbatimTextOutput(session$ns("mergeValue")) )),
                
                fluidRow(
                    
                    column(12,
                           selectInput(inputId = session$ns("scale"),
                                       label = .addBSpopify(label = 'Scale by:', 
                                                            content = paste0("By default for proteomics or metabolomics data,", 
                                                                             "coseq is done onto Z-scores (data scaled by proteins or metabolites) to group them",
                                                                             "according to their expression profile rather than abundance")),
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
                           selectInput(inputId = session$ns("GaussianModel"),
                                       label = .addBSpopify(label = 'Gaussian Model:', 
                                                            content = paste0("For proteomics or metabolomics data, coseq analysis may fail", 
                                                                             "with default GaussianModel parameter.",
                                                                             "In this case, an error message will indicate to switch the other option: Gaussian_pk_Lk_Bk")),
                                       choices  = Gaussian,
                                       selected = Gaussian[1]))
                ),
                fluidRow(
                    column(8,
                           sliderInput(session$ns("K.values"), label = .addBSpopify(label = 'K ranges:', content = "K=number of clusters"), min=2, max=30, value=c(2,7), step=1)),
                    column(4,
                           numericInput(inputId = session$ns("iter"), label=.addBSpopify(label = 'Iteration:', content = "Number of replicates"), value=5, min = 5, max=20, step = 5))
                ),
                hr(),
                fluidRow(
                    
                    column(8,
                           materialSwitch(inputId = session$ns("clustermqCoseq"), label = .addBSpopify(label = 'use remote cluster', content = "send calculation to the cluster"), value = FALSE, status = "success")
                    ),
                    column(4, actionButton(session$ns("runCoSeq"),"Run")))
            ))
    })
    
    #get list of DGE to process
    DEG_list <- reactive({
        getDEList(object = session$userData$FlomicsMultiAssay[[dataset]],
                  contrasts = input$select, operation = input$unionInter)})
    
    # display nbr of selected genes
    output$mergeValue <- renderText({
        
        if(rea.values[[dataset]]$diffValid == FALSE) return()
        
        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
        paste(length(getDEList(object    = dataset.SE,
                        contrasts = input$select, 
                        operation = input$unionInter)), 
              .omicsDic(dataset.SE)$variableName)
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
    
    ##================================ RUN ==================================##
    # run coexpression analysis
    # coseq
    observeEvent(input$runCoSeq, {
        
        # check if no selected DGE list
        if(length(input$select) == 0){
            
            showModal(modalDialog( title = "Error message", 
                                   "Please select at least 1 DEG list."))
        }
        validate({
            need(length(input$select) != 0, 
                 message="Please select at least 1 DEG list")
        })
        
        # check the number of features
        if(length(DEG_list()) < 100){
            
            showModal(modalDialog(title = "Error message", 
                                  paste0("Need at least 100 ", 
                                         .omicsDic(session$userData$FlomicsMultiAssay[[dataset]])$variableName, ".")))
        }
        validate({
            need(length(DEG_list()) >= 100, 
                 message=paste0("Need 100 at least ", 
                                .omicsDic(session$userData$FlomicsMultiAssay[[dataset]])$variableName, "."))
        })
        
        
        # dble execution
        param.list <- list(method.        = "coseq",
                           model          = input$model,
                           gene.list.names= input$select,
                           merge.type     = input$unionInter,
                           K.range        = input$K.values[1]:input$K.values[2],
                           replicates.nb  = input$iter,
                           transformation = input$transfo,
                           normFactors    = input$norm,
                           GaussianModel  = input$GaussianModel, 
                           clustermq      = input$clustermqCoseq,
                           scale          = input$scale)
        
        if(check_run_coseq_execution(
          session$userData$FlomicsMultiAssay[[dataset]], 
          param.list) == FALSE) return()
        
        # initialize reactive value
        rea.values[[dataset]]$coExpAnal  <- FALSE
        rea.values[[dataset]]$coExpAnnot <- FALSE
        rea.values[[dataset]]$CoExpClusterNames <- NULL
        
        # initialize MAE object
        session$userData$FlomicsMultiAssay <-
          resetRflomicsMAE(session$userData$FlomicsMultiAssay,
                           datasetNames = dataset,
                           singleAnalyses = c("CoExpAnal",
                                        "CoExpEnrichAnal"))
        
        #---- progress bar ----#
        progress <- shiny::Progress$new()
        progress$set(message = "Run coseq: ", value = 0)
        on.exit(progress$close())
        progress$inc(1/2, detail = "In progress")
        #----------------------#
        
        message("# 9- CoExpression analysis... ", dataset )
        
        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
        # run coseq
        dataset.SE <- runCoExpression(object = dataset.SE, 
                                      merge = input$unionInter, 
                                      contrastNames = input$select,
                                      K = input$K.values[1]:input$K.values[2], 
                                      replicates = input$iter,
                                      model  = input$model, 
                                      transformation = input$transfo, 
                                      scale = input$scale,
                                      normFactors = input$norm,
                                      GaussianModel = input$GaussianModel, 
                                      clustermq = input$clustermqCoseq, 
                                      cmd = TRUE, silent = TRUE)
        
        session$userData$FlomicsMultiAssay[[dataset]] <- dataset.SE
        
        # If an error occured
        if(isFALSE(metadata(dataset.SE)$CoExpAnal[["results"]])){
            
            showModal(modalDialog( title = "Error message", 
                                   paste0("No results! ", 
                                          as.character(metadata(dataset.SE)$CoExpAnal[["error"]]))))
        }
        
        validate(
            need(!isFALSE(metadata(dataset.SE)$CoExpAnal[["results"]]), 
                 paste0("No results!", as.character(metadata(dataset.SE)$CoExpAnal[["error"]])))
        )
        
        #---- progress bar ----#
        progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
        #----------------------#
        
        rea.values[[dataset]]$coExpAnal  <- TRUE
        rea.values[[dataset]]$CoExpClusterNames <- 
          names(metadata(dataset.SE)$CoExpAnal$clusters)
        
        rea.values$datasetCoEx <-
          getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                  analyses = "CoExpAnal")
        rea.values$datasetCoExAnnot <-
          getAnalyzedDatasetNames(session$userData$FlomicsMultiAssay,
                                  analyses = "CoExpEnrichAnal")
        
    }, ignoreInit = TRUE)
    
    ##============================== DISPLAY ================================##
    
    output$CoExpResultUI <- renderUI({
        
        if (rea.values[[dataset]]$coExpAnal == FALSE) return()
        
        # Extract results
        MAE.data <- session$userData$FlomicsMultiAssay
        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
        
        factors.bio <- getBioFactors(MAE.data)
        
        plot.coseq.res <- plotCoExpression(dataset.SE)
        
        nb_cluster     <- metadata(dataset.SE)$CoExpAnal[["cluster.nb"]]
        coseq.res      <- metadata(dataset.SE)$CoExpAnal[["coseqResults"]]
        cluster.comp   <- metadata(dataset.SE)$CoExpAnal[["clusters"]]
        topDEF         <- metadata(dataset.SE)$DiffExpAnal$mergeDEF
        
        
        box(title = paste0("Number of clusters: ", nb_cluster), 
            status = "warning", solidHeader = TRUE, width = 14,
            
            tabBox( id = "runClustering", width = 12,
                    
                    tabPanel("ICL",  
                             br(),br(),
                             tags$i("Integrated Completed Likelihood (ICL) plotted versus number of clusters. 
                                The ICL is the criterion used to select the number of clusters. 
                                The number of clusters (K) that minmizes the ICL is the best number of clusters."),
                             br(),hr(),br(),
                             renderPlot({ plot.coseq.res$ICL })),
                    tabPanel("probapost_barplots",
                             br(),
                             br(),
                             tags$i(
                                 paste0(
                                     "Barplots giving a quality control of each cluster.",
                                     " In green: number of " ,
                                     .omicsDic(dataset.SE)$variableName,
                                     " that are member of a cluster with a high confidance (Max Conditional Probability  > 0.8).",
                                     " In purple:  number of ",
                                     .omicsDic(dataset.SE)$variableName,
                                     " that are member of a cluster with a low confidance (Max Conditional Probability < 0.8).",
                                     " A high proportion of green observations may indicate that the clustering is not reliable and can't be exploited."
                                 )
                             ),
                             br(),
                             hr(),
                             br(),
                             renderPlot({
                                 plot.coseq.res$probapost_barplots
                             })),
                    tabPanel(
                        "profiles",
                        br(),
                        br(),
                        tags$i(
                            "Overview of cluster's expression profiles. Boxplot are colored acording to  biological  factors"
                        ),
                        br(),
                        hr(),
                        br(),
                        renderPlot({
                            plot.coseq.res$boxplots + 
                                theme(axis.text.x = element_text(angle = 90, hjust = 1))
                        })
                    ),
                    tabPanel("profiles per cluster",
                             fluidRow(
                                 column(3, 
                                        br(),br(),br(),
                                        ## select cluster to plot
                                        box(width = 14, background = "light-blue", title = NULL,
                                            radioButtons(inputId = session$ns("selectCluster"), label = "Select cluster:",
                                                         choices  = seq_len(nb_cluster), selected = 1),
                                            radioButtons(inputId = session$ns("profile.condition"), label = "Condition:",
                                                         choices  = c("groups", factors.bio), selected = "groups"),
                                            uiOutput(session$ns("observationsUI"))
                                        )
                                 ),
                                 column(9, 
                                        br(),br(),
                                        tags$i(tags$p("boxplots of expression/abundance profiles in each cluster colored by experimental factor's modality.")),
                                        tags$i(tags$p("red line join the median of the expression/abundance profiles of each factor's modality")),
                                        br(), hr(),br(),
                                        renderPlot({ 
                                            plotCoExpressionProfile(dataset.SE, input$selectCluster, condition=input$profile.condition, features=input$observations) }))
                             )),
                    
                    tabPanel("cluster composition",
                             fluidRow(
                                 br(),br(),
                                 tags$i(paste0("Cluster's composition according to the ", 
                                               .omicsDic(dataset.SE)$variableName,
                                               "'s contrast belonging")),
                                 br(), hr(),br(),
                                 renderPlot({ 
                                     H <- getCoexpSettings(dataset.SE)$gene.list.names
                                     
                                     if(length(H) > 1){
                                         plotCoseqContrasts(dataset.SE)
                                     }
                                 })
                             ))
            )
        )
    })
    
    output$observationsUI  <- renderUI({
        
        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
        coseq.res  <- metadata(dataset.SE)$CoExpAnal[["coseqResults"]]
        clustr_num <- paste0("Cluster_",input$selectCluster)
        assays.data <- filter(as.data.frame(coseq.res@assays@data[[1]]), 
                              get(clustr_num) > 0.8)
        
        choices <- rownames(assays.data)
        names(choices) <- paste0(choices, " (",assays.data[,clustr_num], ")")
        
        selectizeInput(
            inputId = session$ns("observations"), label = "Observations(prob)",
            choices = choices, multiple = FALSE)
    })
    
    ## summary
    output$ERROR.UI <- renderUI({
        
        dataset.SE <- session$userData$FlomicsMultiAssay[[dataset]]
        
        if(rea.values[[dataset]]$coExpAnal == FALSE || dim(metadata(dataset.SE)$CoExpAnal[["stats"]])[1] == 0) return()
        
        box(title = "Failed cases", width = 14, status = "warning", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
            
            DT::renderDataTable(DT::datatable(as.data.frame(metadata(dataset.SE)$CoExpAnal[["stats"]]), 
                                              options = list(dom = 'tip'), rownames = FALSE))
        )
    })
}

############## functions ###############

# ----- check run diff execution ------
check_run_coseq_execution <- function(object.SE, param.list = NULL){
    
    # filtering setting
    if(length(metadata(object.SE)[["CoExpAnal"]]) == 0 || metadata(object.SE)[["CoExpAnal"]]$results != TRUE) return(TRUE)
    
    if(isFALSE(dplyr::setequal(param.list$gene.list.names, getCoexpSettings(object.SE)$gene.list.names))) return(TRUE)
    if(param.list$model          != getCoexpSettings(object.SE)$model)                       return(TRUE)
    if(param.list$merge.type     != getCoexpSettings(object.SE)$merge.type)                  return(TRUE)
    if(isFALSE(dplyr::setequal(param.list$K.range, getCoexpSettings(object.SE)$K.range)))    return(TRUE)
    if(param.list$replicates.nb  != getCoexpSettings(object.SE)$replicates.nb)               return(TRUE)
    if(param.list$transformation != getCoexpSettings(object.SE)$transformation)              return(TRUE)
    if(param.list$normFactors    != getCoexpSettings(object.SE)$normFactors)                 return(TRUE)
    if(param.list$GaussianModel  != getCoexpSettings(object.SE)$GaussianModel)               return(TRUE)
    if(param.list$scale          != getCoexpSettings(object.SE)$scale)                       return(TRUE)
    
    return(FALSE)
}


