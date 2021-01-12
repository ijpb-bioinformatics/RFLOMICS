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
        box(title = "parameters", solidHeader = TRUE, status = "warning", width = 14,
            
          column(6,
            selectInput(ns("methode"), label = "Method :", 
                        choices = list("normal"="normal"), selected = "normal"),
            selectInput(ns("transfo"), label = "Transformation :", 
                        choices = list("arcsin"="arcsin"), selected = "arcsin"),
            selectInput(ns("norm"), label = "normFactors :", 
                        choices = list("TMM"="TMM"), selected = "TMM")
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
            tabPanel("probapost_boxplots",  plotOutput(ns("probapost_boxplots"))),
            tabPanel("probapost_barplots",  plotOutput(ns("probapost_barplots"))),
            tabPanel("probapost_histogram", plotOutput(ns("probapost_histogram"))),
            tabPanel("boxplots_bis", 
                     uiOutput(ns("selectClusters"))#,
                     #plotOutput(ns("boxplots_bis"))
                     
                     )
          )
        )
      )#,
    # fluidRow(
    #   box(title = "Explore clusters", width = 12, status = "warning")
    # 
    # )
    )
  )
}

CoSeqAnalysis <- function(input, output, session, dataset){
 
  # get coseq parameters
  # 
  # Select lists of DEG to co-expression analysis 
  output$selectDEGtoCoExp <- renderUI({

    checkboxGroupInput(inputId = session$ns("select"), label = "Select DEG to Coseq :",
                       choiceNames  = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                       choiceValues = FlomicsMultiAssay@metadata$design@Contrasts.Sel$tag,
                       selected     = FlomicsMultiAssay@metadata$design@Contrasts.Sel$tag)

  })
  
  # get list of DEG to process 
  # from union or intersection
  
  
  
  #observeEvent(input$select, {
    
  #})
  
  
  #run coseq when button is clicked
  observeEvent(input$runCoSeq, {
   
    progress <- shiny::Progress$new()
    progress$set(message = "Making plot", value = 0)
    on.exit(progress$close())
    
    progress$inc(0, detail = paste("Doing part ", 0,"%", sep=""))
    
    DEG_list <- getDEGlist_for_coseqAnalysis(matrix    = FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg.mat"]], 
                                             colnames  = input$select, 
                                             mergeType = input$unionInter)
    
    output$mergeValue <- renderText({ print(paste(length(DEG_list), "genes", sep =" ")) })
    
    # run coseq
    progress$inc(1/4, detail = paste("Doing part ", 0,"%", sep=""))
 
    # FlomicsMultiAssay <- runCoseq(counts = assay(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]])[DEG_list,] , 
    #                       K=input$minK:input$maxK, iter = input$iter,
    #                       model  = input$methode, transformation=input$transfo, 
    #                       parallel=TRUE, meanFilterCutoff=50, normFactors="TMM")
    
    FlomicsMultiAssay <<- runCoExpression(FlomicsMultiAssay, data = paste0(dataset,".filtred"), "coseq", DEG_list, 
                                  K=input$minK:input$maxK, iter = input$iter, model  = input$methode, 
                                  transformation=input$transfo, normFactors="TMM")
    
    
    
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    
    coseq.res <- FlomicsMultiAssay@ExperimentList[[ paste0(dataset,".filtred")]]@metadata[["CoExpResults"]][["coseqResults"]]
    groups    <- unite(as.data.frame(FlomicsMultiAssay@colData[FlomicsMultiAssay@metadata$design@Factors.Type == "Bio"]), 
                       col="groups", sep="_", remove = TRUE) %>% dplyr::mutate(samples = rownames(.)) %>%
                 dplyr::arrange(factor(samples, levels = names(coseq.res@y_profiles)))
    plot.coseq.res <- coseq::plot(coseq.res, conds = groups$groups)
    
    output$logLike  <- renderPlot({ plot.coseq.res$logLike }) 
    output$ICL      <- renderPlot({ plot.coseq.res$ICL }) 
    output$profiles <- renderPlot({ plot.coseq.res$profiles }) 
    output$boxplots <- renderPlot({ plot.coseq.res$boxplots }) 
    output$probapost_boxplots  <- renderPlot({ plot.coseq.res$probapost_boxplots }) 
    output$probapost_barplots  <- renderPlot({ plot.coseq.res$probapost_barplots }) 
    output$probapost_histogram <- renderPlot({ plot.coseq.res$probapost_histogram }) 
    
    
    output$selectClusters <- renderUI({
      
      nb_cluster <- coseq.res@metadata$nbCluster[min(coseq.res@metadata$ICL) == coseq.res@metadata$ICL]
      
      fluidRow(
        
        ## plot selected cluster(s)
        renderPlot({ coseq.y_profile.one.plot(coseq.res, input$selectCluster, groups) }),
        
        ## select cluster to plot
        checkboxGroupInput(inputId = session$ns("selectCluster"), label = "Select cluster(s) :",
                           choices  = 1:nb_cluster, selected = 1:nb_cluster, inline = TRUE)
      )
    })
    
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    
  })
}
