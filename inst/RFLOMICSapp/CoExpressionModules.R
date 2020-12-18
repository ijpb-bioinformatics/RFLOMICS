CoSeqAnalysisUI <- function(id){
  
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    
    ### parametres for Co-Exp
    fluidRow( 
      #uiOutput(ns("DiffParam"))
      box(title = "", width = 12, status = "warning",
          column(6,
                 uiOutput(ns("selectDEGtoCoExp"))
          ),
          column(4,
                 uiOutput(ns("selectMergeOrInters"))
          ),
          column(2,
                 fluidRow( 
                   selectInput(ns("tools"), label = "Tools :",
                               choices = list("coseq"="coseq"),
                               selected = "coseq")
                 ),
                 fluidRow( 
                   selectInput(ns("methode"), label = "Method :",
                               choices = list("normal"="normal"),
                               selected = "normal")
                 ),
                 fluidRow( 
                   selectInput(inputId = ns("clustermq-coseq"),
                               label="send job to cluster",
                               choices = list("no"=FALSE,"genotoul"=TRUE))
                 ),
                 fluidRow( actionButton(ns("runCoSeq"),"Run Co-Expression Analysis") )
          )
      )
      
    ),
    fluidRow( 

    ),
    fluidRow(
      box(title = "Selected 1 :", width = 12, status = "warning",
        verbatimTextOutput(ns("res_select"))
      ),
      box(title = "Selected 2 :", width = 12, status = "warning",
          verbatimTextOutput(ns("unionInter_select"))
      )
      
      
      
    )#,
    # fluidRow(
    #   box(title = "Selected 1 :", width = 12, status = "warning",
    #       plotOutput(ns("upset"))
    #   )
    # )
  )
}

CoSeqAnalysis <- function(input, output, session, dataset){
 
  # Select lists of DEG to co-expression analysis 
  output$selectDEGtoCoExp <- renderUI({
  
    checkboxGroupInput(inputId = session$ns("select"), label = "Select DEG to Coseq :", 
                       choiceNames  = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
                       choiceValues = FlomicsMultiAssay@metadata$design@Contrasts.Sel$tag,
                       selected     = FlomicsMultiAssay@metadata$design@Contrasts.Sel$tag)
  })
  
  # Select type of merge : union or intersection
  output$selectMergeOrInters <- renderUI({
    
    
    
    # checkboxGroupInput(inputId      = session$ns("unionInter"), label = "union or intersection", 
    #                    choiceNames  = c(paste("union",        dim(dplyr::filter(merge_DEG, rowSums != 0))[1], sep=""), 
    #                                     paste("intersection", dim(dplyr::filter(merge_DEG, rowSums == 3))[1], sep="")),
    #                    choiceValues = c("union", "intersection"),
    #                    selected    = "union")
    
    radioButtons(inputId      = session$ns("unionInter"), label = "union or intersection", 
                 choiceNames  = c("union", "intersection"),
                 choiceValues = c("union", "intersection"),
                 selected    = "union")
  })
  
  output$res_select <- renderPrint({
    input$select
  })
  
  output$unionInter_select <- renderPrint({
    
    merge_DEG <- dplyr::select(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg.mat"]], input$select) %>%
      dplyr::mutate(rowSums = rowSums(.))
    
      switch(input$unionInter ,
             "union"={ dim(dplyr::filter(merge_DEG, rowSums != 0))[1] },
             "intersection"={ dim(dplyr::filter(merge_DEG, rowSums == length(input$select)))[1] }
      )
  })
  
  # union_DEG <- dplyr::select(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg.mat"]], select) %>%
  #              dplyr::mutate(rowSums = rowSums())
  
  # output$upset <- renderUI({
  # 
  #     renderPlot({
  #       DEG_mat <- dplyr::select(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg.mat"]], select)
  #       UpSetR::upset(DEG_mat, sets = (names(DEG_mat[,-1])))
  #     })
  # })
  
  
}
