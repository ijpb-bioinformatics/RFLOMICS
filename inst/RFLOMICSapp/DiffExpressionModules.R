


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
                     selectInput(ns("AnaDiffMethod"), label = "Method :",
                                 choices = list("glmfit (edgeR)"="edgeRglmfit"),
                                 selected = "edgeRglmfit")
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

    progress <- shiny::Progress$new()
    progress$set(message = "Run Diff", value = 0)
    on.exit(progress$close())
    
    # run diff analysis with select method
    progress$inc(1/10, detail = paste("Doing part ", 10,"%", sep=""))
    FlomicsMultiAssay <<- RunDiffAnalysis(FlomicsMultiAssay, data=paste0(dataset,".filtred"), 
                                          contrastList = input$contrastList, FDR =input$FDRSeuil, 
                                          DiffAnalysisMethod=input$AnaDiffMethod,
                                          clustermq=input$clustermq)
    
    progress$inc(1/2, detail = paste("Doing part ", 50,"%", sep=""))
    
    Contrasts.Sel <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["contrasts"]]
    
    output$ContrastsResults <- renderUI({
      
      lapply(1:length(Contrasts.Sel$contrast), function(i) {
        
        vect     <- unlist(Contrasts.Sel[i,])
        res      <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["DGELRT"]][[vect["contrastName"]]]
        resTable <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["TopDGE"]][[vect["contrastName"]]]

        fluidRow(
          column(10,
                 box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = paste0(vect["tag"], " : ", vect["contrastName"], sep=""),

                     
                     column(6,
                            ### MAplot
                         renderPlot({
                           res.FDR <-topTags(res, n = dim(res)[1])
                           MA.plot(data = res.FDR$table, FDRcutoff = input$FDRSeuil,
                                   pngFile =file.path(tempdir(), paste0(dataset, "_MAplot_", gsub(" ", "", vect["contrastName"]), ".png")))
                         })
                     ),
                     column(6,
                            ### pvalue plot ###
                            renderPlot({
                              
                              pvalue.plot(data    =resTable,
                                          contrast=vect["contrastName"],
                                          pngFile =file.path(tempdir(), paste0(dataset, "_PvalueDistribution_", gsub(" ", "", vect["contrastName"]), ".png")))
                            })
                     ),
                     verticalLayout(
                       
                       
                       tags$br(),
                       ### DEG result table ###
                       DT::renderDataTable({
                         
                         DT::datatable(round(resTable[resTable$FDR <= input$FDRSeuil,],5),
                                       options = list(rownames = FALSE, pageLength = 10))
                       })
                       
                     )
                 )
          )
        )
      })
    })
    progress$inc(3/4, detail = paste("Doing part ", 75,"%", sep=""))
    

    ### intersection
    DEG_mat <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["mergeDGE"]]

      output$ResultsMerge <- renderUI({
              if (ncol(DEG_mat) > 2){ 
                    fluidRow(
                      column(10,
                             box(width=12,  status = "warning", 
                                 
                                 renderPlot({ UpSetR::upset(DEG_mat, sets = (names(DEG_mat[,-1]))) })
                                 )
                             )
                      )
              }else{
                fluidRow(
                  
                )
              }
        })
    

    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
  })
  
  return(input)
}


