AnnotationEnrichmentUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    fluidRow(
      box(title = span(tagList(icon("cogs"), "   Enrichment analysis")),
          solidHeader = TRUE, status = "warning", width = 12, )
      ),
    ## parameters + input data
    fluidRow(
      
      box(title = span(tagList(icon("cogs"), "")), solidHeader = FALSE, status = "warning", width = 3, 

        ## gene lists
        uiOutput(ns("selectGeneListtoAnnot")),
        
        ## file with annotation (geneID, Term, Name, Domain/source)
        fileInput(inputId = ns("annotationFile"), label = "Annotation file :", 
                            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        
        ## choice of probability
        selectInput(ns("EnrichMethod"), label = "Test :", 
                           choices = list("Hypergeometric"="hypergeometric"),selected = "hypergeometric"),
        
        ## alpha threshold
        numericInput(inputId = ns("Alpha_Enrichment"), label="P-value :", value=0.01 , min = 0, max=0.1, step = 0.01),
                  actionButton(ns("runEnrich"),"Run")
        
      ),
      uiOutput(ns("AnnotEnrichResults"))
    )

  )
}

AnnotationEnrichment <- function(input, output, session, dataset){
  
  output$selectGeneListtoAnnot <- renderUI({
    
    ListNames.diff <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["contrasts"]]$contrastName
    ListNames.coseq <- names(FlomicsMultiAssay@ExperimentList[[ paste0(dataset,".filtred")]]@metadata$CoExpAnal[["clusters"]])
    
    # multiInput(
    #   inputId = session$ns("GeneList.diff"),
    #   label = "Select DEG to Coseq :", 
    #   choiceNames = ListNames,
    #   choiceValues = ListNames
    # )
    list(
      pickerInput(
        inputId = session$ns("GeneList.diff"),
        label = "Select DEG lists:", 
        choices = ListNames.diff,
        options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
        multiple = TRUE, selected = ListNames.diff
      ),
      
      pickerInput(
        inputId = session$ns("GeneList.coseq"), label = "Select Clusters :",
        choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
        options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        
      )
    )
  })
  
  
  
  
  observeEvent(input$runEnrich, {
    
    # check list of genes
    if(length(c(input$GeneList.diff, input$GeneList.coseq )) == 0){
      
      showModal(modalDialog( title = "Error message", "Please select at least 1 gene list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff, input$GeneList.coseq )) != 0, message="Please select at least 1 gene list") 
    })
    
    # check annotation file
    if(is.null(input$annotationFile$datapath)){
      
      showModal(modalDialog( title = "Error message", "need annotation file."))
    }
    validate({ 
      need(!is.null(input$annotationFile$datapath), message="need annotation file") 
    })
    
    annotation <- fread(file = input$annotationFile$datapath, sep="\t", header = TRUE)
    colnames(annotation) <- c("geneID", "Term", "Name", "Domain")
    
    
    print(paste("# 11- Enrichment Analysis...", dataset))
    
    ## run annotation
    dataset.SE <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]
    dataset.SE <- runAnnotationEnrichment(dataset.SE, annotation= annotation, 
                                          DiffListNames=input$GeneList.diff, CoExpListNames=input$GeneList.coseq, 
                                          alpha = input$Alpha_Enrichment, probaMethod = input$EnrichMethod)
    FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]] <<- dataset.SE
    
    ## print results
    output$AnnotEnrichResults <- renderUI({

      # foreach gene list selected (contrast)
      lapply(names(dataset.SE@metadata$EnrichAnal[["results"]]), function(listname) {

        data <- dataset.SE@metadata$EnrichAnal[["results"]][[listname]][["Over_Under_Results"]]
        
        fluidRow(
          column(9, 
                 box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname,
                     
                     ### pvalue plots
                     tabBox( id = "pvaluePlot", width = 12,
                             
                             tabPanel(title = "overrepresented", 
                                      
                                  column(10,
                                        
                                        renderPlot({ Enrichment.plot(dataset.SE, Over_Under="overrepresented", top = input$top.over, listNames=listname) })
                                  ),
                                  column(2,
                                         
                                        radioGroupButtons(direction = "vertical", inputId = session$ns("top.over"), 
                                                          label = "", choices = c("all",  50), justified = FALSE, selected = 50)
                                  ),
                                  verticalLayout(
                                    
                                    tags$br(),
                                    ### DEG result table ###
                                    DT::renderDataTable({
                                      
                                      DT::datatable(dplyr::filter(data, Decision == "overrepresented"), 
                                                    options = list(rownames = FALSE, pageLength = 5, lengthMenu = c(5,10, 15, 20), scrollX = T))
                                    })
                                  )
                                      
                             ),
                             tabPanel(title = "underrepresented", 
                                      
                                  column(10,
                                         renderPlot({Enrichment.plot(dataset.SE, Over_Under="underrepresented", top = input$top.under, listNames=listname)})
                                  ),
                                  column(2,
                                         radioGroupButtons(direction = "vertical", inputId = session$ns("top.under"), 
                                                           label = "", choices = c("all",  50), justified = FALSE , selected = 50)
                                  ),
                                  verticalLayout(
                                    
                                    tags$br(),
                                    ### DEG result table ###
                                    DT::renderDataTable({
                                      
                                      DT::datatable(dplyr::filter(data, Decision == "underrepresented"), 
                                                    options = list(rownames = FALSE, pageLength = 5, lengthMenu = c(5,10, 15, 20), scrollX = T))
                                    })
                                  )
                             )
                        )
                     )
                 )
          )
        })
      })
    })

}
    
    