AnnotationEnrichmentUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    
    ## parameters + input data
    fluidRow(
      
      box(title = span(tagList(icon("cogs"), "")), solidHeader = FALSE, status = "warning", width = 12, 

        ## gene lists
        column(4, uiOutput(ns("selectGeneListtoAnnot"))),
        
        ## file with annotation (geneID, Term, Name, Domain/source)
        column(3, fileInput(inputId = ns("annotationFile"), label = "Annotation file :", 
                            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
        
        ## choice of probability
        column(2, 
               selectInput(ns("EnrichMethod"), label = "Test :", 
                           choices = list("Hypergeometric"="hypergeometric"),selected = "hypergeometric")),
        
        ## alpha threshold
        column(2, numericInput(inputId = ns("Alpha_Enrichment"), label="P-value :", value=0.01 , min = 0, max=0.1, step = 0.01),
                  actionButton(ns("runEnrich"),"Run"))
        
      )
    ),
        
    fluidRow( uiOutput(ns("AnnotEnrichResults")) )
    
  )
}

AnnotationEnrichment <- function(input, output, session, dataset){
  
  output$selectGeneListtoAnnot <- renderUI({
    
    ListNames.diff <- FlomicsMultiAssay@metadata$design@Contrasts.Sel$contrastName
    ListNames.coseq <- names(FlomicsMultiAssay@ExperimentList[[ paste0(dataset,".filtred")]]@metadata[["CoExpResults"]][["clusters"]])
    
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
    
    
    print("# 11- Enrichment Analysis...")
    
    ## list of gene list to annotate
    geneLists <- list()
    geneLists.diff <- list()
    geneLists.diff <- lapply(input$GeneList.diff, function(listname){

       row.names(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]@metadata$AnaDiffDeg[[listname]])
      
    })
    names(geneLists.diff) <- input$GeneList.diff
    
    geneLists.coseq <- list()
    geneLists.coseq <- lapply(input$GeneList.coseq, function(listname){
      
      FlomicsMultiAssay@ExperimentList[[ paste0(dataset,".filtred")]]@metadata[["CoExpResults"]][["clusters"]][[listname]]
      
    })
    names(geneLists.coseq) <- input$GeneList.coseq
    
    geneLists <- c(geneLists.diff, geneLists.coseq)
    
    ## run annotation
    FlomicsMultiAssay <<- runAnnotationEnrichment(FlomicsMultiAssay, data = paste0(dataset,".filtred"), 
                                                  annotation, geneLists, alpha = input$Alpha_Enrichment, probaMethod = input$EnrichMethod)

    ## print results
    output$AnnotEnrichResults <- renderUI({

      # foreach gene list selected (contrast)
      lapply(names(FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]$metadata$AnnotEnrich), function(listname) {

        data <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]$metadata$AnnotEnrich[[listname]][["Over_Under_Results"]]
        
        fluidRow(
          column(12, 
                 box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname,
                     
                     ### pvalue plots
                     tabBox( id = "pvaluePlot", width = 12,
                             
                             tabPanel(title = "overrepresented", 
                                      
                                  column(10,
                                        
                                        renderPlot({ pvalue.enrichment.plot(data, "overrepresented", index = input$top.over, pngFile=NULL)})
                                  ),
                                  column(2,
                                         
                                        radioGroupButtons(direction = "vertical", inputId = session$ns("top.over"), 
                                                          label = "", choices = c("all",  "Top50"), justified = FALSE, selected = "Top50")
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
                                         renderPlot({ pvalue.enrichment.plot(data, "underrepresented", index = input$top.under, pngFile=NULL)})
                                  ),
                                  column(2,
                                         radioGroupButtons(direction = "vertical", inputId = session$ns("top.under"), 
                                                           label = "", choices = c("all",  "Top50"), justified = FALSE , selected = "Top50")
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
    
    