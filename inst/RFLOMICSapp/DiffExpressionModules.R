


DiffExpAnalysisUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    
    ### parametres for Diff Analysis
    fluidRow( 
      #uiOutput(ns("DiffParam"))
      box(title = "", width = 12, status = "warning",
          column(5,
                 selectInput(ns("AnaDiffMethod"), label = "Method :",
                             choices = list("glmfit (edgeR)"="edgeRglmfit",
                                            "limma" = "limma"),
                             selected = "edgeRglmfit")
          ),
          column(3,
                 numericInput(inputId = ns("FDRSeuil"),
                              label="FDR :",
                              value=0.05, 0, max=1, 0.01)
          ),
          column(3,
                 selectInput(inputId = ns("clustermq"),
                             label="send job to cluster",
                             choices = list("no"=FALSE,"genotoul"=TRUE))
          ),
          column(5,
                 actionButton(ns("runAnaDiff"),"Run the differential analysis")
          )
      )
      
    ),
    tags$br(),
    tags$br(),
    fluidRow( 
      
      uiOutput(ns("ContrastsResults"))
      
      
      )#,
    #fluidRow( uiOutput(ns("ResultsMerge")))
  )
}



DiffExpAnalysis <- function(input, output, session, dataset){
  
 
  # Run the differential analysis for each contrast set
  #   -> return a dynamic user interface with a collapsible box for each contrast
  #         - Pvalue graph
  #         - MAplot
  #         - Table of the DE genes 
  observeEvent(input$runAnaDiff, {
    
    print("# 9- Diff Analysis...")

    # run diff analysis with select method
    FlomicsMultiAssay <<- RunDiffAnalysis(FlomicsMultiAssay, data=paste0(dataset,".filtred"),
                                          FDR =input$FDRSeuil , DiffAnalysisMethod=input$AnaDiffMethod,
                                          clustermq=input$clustermq)

    output$ContrastsResults <- renderUI({

      vect <- as.vector(FlomicsMultiAssay@metadata$design@Contrasts.List$hypoth)
      names(vect) <- as.vector(FlomicsMultiAssay@metadata$design@Contrasts.List$idContrast)

      lapply(FlomicsMultiAssay@metadata$design@Contrasts.Sel, function(i) {

        resTable <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg"]][[i]]



        fluidRow(
          column(10,
                 box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = paste0(i, " : ", vect[i]),

                     verticalLayout(

                       ### pvalue plot ###
                       renderPlot({

                         pvalue.plot(data    =FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg"]][[i]],
                                     contrast=as.vector(filter(FlomicsMultiAssay@metadata$design@Contrasts.List , idContrast == i)$hypoth),
                                     pngFile =file.path(tempdir(), paste0(dataset, "_PvalueDistribution_", gsub(" ", "", vect[i]), ".png")))
                       }),
                       tags$br(),
                       DT::renderDataTable({

                         DT::datatable(round(resTable[resTable$FDR <= input$FDRSeuil,],5),
                                       options = list(rownames = FALSE, pageLength = 10))
                       })

                     )
                 )
          ),
          column(2, checkboxInput(inputId = "checkContrasts", label = "validate" ,value = TRUE , width=1))
        )
      })
    })
    
  })
  
}



# DiffExpMergeUI <- function(id){
#   
#   #name space for id
#   ns <- NS(id)
#   
#   tagList(  
#     fluidRow( uiOutput(ns("ResultsMerge")))
#   )
# }
# 
# 
# DiffExpMerge <- function(input, output, session, dataset){
#   
#   # merge diff results
#   mat2venn <- list()
#   for(i in FlomicsMultiAssay@metadata$design@Contrasts.Sel) {
#     
#     mat2venn[[i]][["features"]] <-  row.names(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg"]][[i]])
#     mat2venn[[i]][[i]] <- rep(1, dim(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata[["AnaDiffDeg"]][[i]])[1])
#     mat2venn[[i]] <- tbl_df(mat2venn[[i]])
#   }
#   
#   mat2venn.df <- mat2venn %>% purrr::reduce(dplyr::full_join, by="features")
#   
#   mat2venn.df[is.na(mat2venn.df)] <- 0
#   
#   output$ResultsMerge <- renderUI({
#     fluidRow(
#       column(10,
#              box(width=12, solidHeader = TRUE, title = "Combine results", collapsible = TRUE, collapsed = TRUE, status = "warning",
#                  column(width = 7,
#                         renderPlot({
#                           title(main = "snp")
#                           ven::venn(mat2venn.df[,-1] , ilab=TRUE, zcolor = "style")
#                         })
#                  ),
#                  column(width = 5,
#                         radioButtons("choise", label="" , choices = c("union","intersect"), selected = "union",
#                                      inline = FALSE, width = 2, choiceNames = NULL, choiceValues = NULL),
#                         actionButton("buttonValidMerge","Valid")
#                  )
#              )
#       )
#     )
#   })
#   
#   
# }