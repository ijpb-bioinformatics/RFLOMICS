


DiffExpAnalysisUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    
    ### parametres for Diff Analysis
    fluidRow( 
      #uiOutput(ns("DiffParam"))
      box(title = span(tagList(icon("cogs"), "   edgeR")), width = 12, status = "warning",
          
              column(4,  
                     uiOutput(ns("contrastListUI")),
                     actionButton(ns("runAnaDiff"),"Run")
                     ),
              column(3,
                     uiOutput(ns("AnaDiffMethodUI")),
                     
                     ),
              column(2,
                     numericInput(inputId = ns("FDRSeuil"),
                                  label="FDR :",
                                  value=0.05, 0, max=1, 0.01)
                     ),
              column(3,
                     selectInput(inputId = ns("clustermq"),
                                 label="send job to cluster",
                                 choices = list("no"=FALSE,"genotoul"=TRUE))
                     )
          )
      
    ),
    tags$br(),
    tags$br(),
    fluidRow( uiOutput(ns("ContrastsResults"))),
    fluidRow( uiOutput(ns("ResultsMerge")))
  )
}



DiffExpAnalysis <- function(input, output, session, dataset){
  
  # list of selected contrast
  output$contrastListUI <- renderUI({
      pickerInput(
        inputId = session$ns("contrastList"),
        label = "Selected contrast :", 
        choices = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName,
        options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
        multiple = TRUE, selected = FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName )
  })
  
  output$AnaDiffMethodUI <- renderUI({
  
        MethodList <- c("glmfit (edgeR)"="edgeRglmfit", "Limma"="limma")
        
        method <- switch (FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$omicType,
                  "RNAseq"       = MethodList[1],
                  "proteomics"   = MethodList[2],
                  "metabolomics" = MethodList[2])
        
        selectInput(session$ns("AnaDiffMethod"), label = "Method :",
                    choices  = method,
                    selected = method)
  })
 
  # Run the differential analysis for each contrast set
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
    
    #---- progress bar ----#    
    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = paste("Doing part ", 10,"%", sep=""))
    #----------------------#

    # run diff analysis with select method
    dataset.SE <- RunDiffAnalysis(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]], design = Design,
                                  contrastList = input$contrastList, FDR =input$FDRSeuil, DiffAnalysisMethod=input$AnaDiffMethod, 
                                  clustermq=input$clustermq)
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- dataset.SE
    
    #---- progress bar ----#    
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#
    
    output$ContrastsResults <- renderUI({

      Contrasts.Sel <- dataset.SE@metadata$DiffExpAnal[["contrasts"]]
      
      lapply(1:length(Contrasts.Sel$contrast), function(i) {
        
        vect     <- unlist(Contrasts.Sel[i,])
        res      <- dataset.SE@metadata$DiffExpAnal[["DGELRT"]][[vect["contrastName"]]]
        
        diff.plots <- DiffAnal.plot(dataset.SE, hypothesis=vect["contrastName"])
        
        fluidRow(
          column(10,
             box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = paste0(vect["tag"], " : ", vect["contrastName"], sep=""),
                 
               ### MAplot
               column(6, renderPlot({ diff.plots$MA.plot }) ),
               
               ### pvalue plot ###
               column(6, renderPlot({ diff.plots$Pvalue.hist }) ),
               
               ### DEG result table ###
               verticalLayout(
                 
                 tags$br(),
                 ### DEG result table ###
                 DT::renderDataTable({ 
                   resTable <- dataset.SE@metadata$DiffExpAnal[["TopDGE"]][[vect["contrastName"]]]
                   DT::datatable(round(resTable[resTable$FDR <= input$FDRSeuil,],5), options = list(rownames = FALSE, pageLength = 10)) 
                })
               )
             )
          )
        )
      })
    })
    
    #---- progress bar ----#    
    progress$inc(3/4, detail = paste("Doing part ", 75,"%", sep=""))
    #----------------------#

    ### intersection
    DEG_mat <- dataset.SE@metadata$DiffExpAnal[["mergeDGE"]]

    output$ResultsMerge <- renderUI({
      if (ncol(DEG_mat) > 2){ 
          fluidRow(
            column(10,
               box(width=12,  status = "warning", 
                   
                   renderPlot({ UpSetR::upset(DEG_mat, sets = (names(DEG_mat[,-1]))) })
                   )
               )
          )
        }
      })
    
      #---- progress bar ----#    
      progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
      #----------------------#
  })
  
  return(input)
}


