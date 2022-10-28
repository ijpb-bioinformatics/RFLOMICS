AnnotationEnrichmentClusterProfUI <- function(id){
  
  options(shiny.maxRequestSize = 3000*1024^2)
  
  #name space for id
  ns <- NS(id)
  
  tagList(  
    fluidRow(
      box(title = span(tagList(icon("chart-pie"), "   Enrichment analysis ", tags$small("(Scroll down for instructions)"))),
          solidHeader = TRUE, status = "warning", width = 12, collapsible = TRUE, collapsed = TRUE,
          p("..."))
    ),
    ## parameters + input data
    fluidRow(
      column(3, uiOutput(ns("AnnotParamCPRUI")),
             uiOutput(ns("AnnotParamCPRUI2"))
      ),
      column(9, uiOutput(ns("AnnotDiffResultsCPR"))),
      # column(9, uiOutput(ns("AnnotCoExpResults")))
    )
    
  )
}

AnnotationEnrichmentClusterProf <- function(input, output, session, dataset, rea.values){
  
  local.rea.values <- reactiveValues(dataset.SE = NULL)
  local.rea.values$settings_ok <- NULL
  
  # ---- First set of settings: ----
  output$AnnotParamCPRUI <- renderUI({
    
    validate(
      need(rea.values[[dataset]]$diffValid != FALSE, "Please run diff analysis and validate your choices...")
    )
    
    #### TODO DELETE
    # load("inst/ExamplesFiles/FlomicsMultiAssay.RData")
    # session <- list()
    # session$userData <- list()
    # session$userData$FlomicsMultiAssay <- FlomicsMultiAssay
    # dataset <- c("proteomics.set1")
    #### TODO DELETE
    
    ## gene lists
    ListNames.diff  <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$DiffExpAnal[["Validcontrasts"]]$contrastName
    
    box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
        
        # select DEG list
        pickerInput(
          inputId = session$ns("GeneList.diff"), label = "Select DEG lists:", 
          choices = ListNames.diff, multiple = TRUE, selected = ListNames.diff,
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        # select clusters
        uiOutput(session$ns("GeneList.coseqCPRUI")),
        
        # select enrichDomain
        pickerInput(
          inputId = session$ns("dom.select"), label = "Select Domain:",
          choices = c("custom", "GO", "GO:BP", "GO:CC", "GO:MF", "KEGG"),
          selected = "",
          options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
        ),
        
        ## alpha threshold
        numericInput(inputId = session$ns("Alpha_Enrichment"), label="P-value :", value = 0.01 , min = 0, max = 0.1, step = 0.01),
        
        actionButton(inputId = session$ns("settings_ok"), "set"),
    )
    
  })
  
  
  # ---- Second set of settings: ----
  local.rea.values$settings_ok <- NULL
  observeEvent(input$settings_ok, {local.rea.values$settings_ok <- TRUE})
  
  output$AnnotParamCPRUI2 <- renderUI({
    if(is.null(local.rea.values$settings_ok)) return()
    
    if(input$dom.select == "custom"){
      box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
          fileInput(inputId = session$ns("annotationFile"), label = "Annotation file :", 
                    accept  = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
          actionButton(inputId = session$ns("runEnrich"),"Run")
      )
    }else if(input$dom.select == "KEGG"){
      box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
          # select keytype
          pickerInput(
            inputId = session$ns("keytype"), label = "Select KeyType:",
            choices = c("TAIR", "kegg", "ENTREZID", "SYMBOL"),
            selected = "kegg",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")),
          
          # type organism
          textInput(inputId = session$ns("KEGG_org"), 
                    label = "Organism", value = "ath", 
                    width = NULL, placeholder = NULL),
          
          actionButton(inputId = session$ns("runEnrich"),"Run")
      )
    }else if(length(grep("GO", input$dom.select))!=0){
      box(title = span(tagList(icon("sliders-h"), "  ", "Setting")), width = 14, status = "warning",
          # select database
          pickerInput(
            inputId = session$ns("db.select"), label = "Select Data Base:",
            choices = c("custom", dir(.libPaths())[grep("org.*db", dir(.libPaths()))]),
            selected = "custom",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          
          # select keytype
          pickerInput(
            inputId = session$ns("keytype"), label = "Select KeyType:",
            choices = c("TAIR", "kegg", "ENTREZID", "SYMBOL"),
            selected = "ENTEZID",
            options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
          ),
          actionButton(inputId = session$ns("runEnrich"),"Run")
      )
    }
  })
  
  output$GeneList.coseqCPRUI <- renderUI({

    if(rea.values[[dataset]]$coExpAnal == FALSE) return()

    ListNames.coseq <- names(session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$CoExpAnal[["clusters"]])

    pickerInput(
      inputId = session$ns("GeneList.coseq"), label = "Select Clusters :",
      choices = ListNames.coseq, multiple = TRUE, selected = ListNames.coseq,
      options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3")
    )
  })
  
  
  observeEvent(input$runEnrich, {
    
    #---- progress bar ----#
    progress <- shiny::Progress$new()
    progress$set(message = "Run Annot", value = 0)
    on.exit(progress$close())
    progress$inc(1/10, detail = "in progress...")
    #----------------------#
    
    local.rea.values$dataset.SE <- session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]]
    
    rea.values[[dataset]]$diffAnnot  <- FALSE
    rea.values[[dataset]]$coExpAnnot <- FALSE
    
    local.rea.values$list_arg <- list()
    local.rea.values$func_to_use <- NULL
    
    # check list of genes
    if(length(c(input$GeneList.diff, input$GeneList.coseq )) == 0){
      showModal(modalDialog( title = "Error message", "Please select at least 1 gene list."))
    }
    validate({ 
      need(length(c(input$GeneList.diff, input$GeneList.coseq )) != 0, message="Please select at least 1 gene list") 
    })
    
    # check annotation file
    if(input$db.select == "custom"){
      local.rea.values$func_to_use <- "enricher"
      
      if(is.null(input$annotationFile$datapath)){
        showModal(modalDialog( title = "Error message", "need annotation file."))
      }
      validate({
        need(!is.null(input$annotationFile$datapath), message="need annotation file")
      })
      annotation <- fread(file = input$annotationFile$datapath, sep="\t", header = TRUE)
      
      if(length(intersect(c("geneID", "Term", "Name", "Domain"), colnames(annotation))) < 4){
        
        showModal(modalDialog( title = "Error message", "The header of the annotation file does not match (see input vignette)..."))
      }
      validate({ 
        need(length(intersect(c("geneID", "Term", "Name", "Domain"), colnames(annotation))) == 4, message="The header of the annotation file does not match (see input vignette)...") 
      })
      
      # reduce ref to genes present in matrix count after filtering
      annotation <- dplyr::filter(annotation , geneID %in% names(local.rea.values$dataset.SE))
      
      #check if ref correspond to features in lists
      if(dim(annotation)[1] == 0){
        showModal(modalDialog( title = "Error message", "There is no correspondence between the feature lists and the annotation file"))
      }
      validate({ 
        need(dim(annotation)[1] != 0, message="There is no correspondence between the feature lists and the annotation file") 
      })
      
      
    }else{
      library(input$db.select)
      
      local.rea.values$domain <- input$dom.select
      
      if(length(grep('GO:', input$dom.select))!=0){
        local.rea.values$domain <- str_split(input$dom.select, ":")[[1]][1]
        local.rea.values$list_arg$ont <- str_split(input$dom.select, ":")[[1]][2]
      }
      
      local.rea.values$func_to_use <- paste0("enrich", local.rea.values$domain)
      
      local.rea.values$list_arg$universe <- names(local.rea.values$dataset.SE)
      local.rea.values$list_arg$keyType <- input$keytype
      local.rea.values$list_arg$qvalueCutoff <- input$Alpha_Enrichment
      
    }
    
    print(paste("# 11- Enrichment Analysis...", dataset))
    
    #---- progress bar ----#
    progress$inc(1/5, detail = paste("Doing part ", 50,"%", sep=""))
    #----------------------#
    
    ## run annotation
    
    if(length(input$GeneList.diff) != 0){
      
      # local.rea.values$dataset.SE <- runAnnotationEnrichment(local.rea.values$dataset.SE, annotation= annotation, 
      #                                                        ListNames=input$GeneList.diff, from="DiffExpAnal", 
      #                                                        alpha = input$Alpha_Enrichment, probaMethod = input$EnrichMethod)
      
      
      
      local.rea.values$resAnnot <- NULL 
      
      rea.values[[dataset]]$diffAnnot <- TRUE
    }
    
    if(length(input$GeneList.coseq) != 0){
      
      local.rea.values$dataset.SE <- runAnnotationEnrichment(local.rea.values$dataset.SE, annotation= annotation, 
                                                             ListNames=input$GeneList.coseq, from="CoExpAnal", 
                                                             alpha = input$Alpha_Enrichment, probaMethod = input$EnrichMethod)
      
      rea.values[[dataset]]$coExpAnnot <- TRUE
    }
    
    session$userData$FlomicsMultiAssay[[paste0(dataset,".filtred")]] <- local.rea.values$dataset.SE
    
    
    #---- progress bar ----#
    progress$inc(1, detail = paste("Doing part ", 100,"%", sep=""))
    #----------------------#
  })
  
  ## print results
  output$AnnotDiffResultsCPR <- renderUI({
    
    if(rea.values[[dataset]]$diffAnnot == FALSE) return()
    
    #dataset.SE <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]
    
    # check if result is empty
    if(is.null(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]])){
      
      showModal(modalDialog( title = "Error message", "There is no correspondence between the feature lists and the annotation file"))
    }
    validate({ 
      need(!is.null(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]]), message="There is no correspondence between the feature lists and the annotation file") 
    })
    
    # foreach gene list selected (contrast)
    lapply(names(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]]), function(listname) {
      
      # if result is empty
      if(is.null(local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]][[listname]])){
        
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
            "There is no correspondence between the feature list and the annotation file !"
        )
        
      }else{
        
        # display result for each list
        data <- local.rea.values$dataset.SE@metadata$DiffExpEnrichAnal[["results"]][[listname]][["Over_Under_Results"]]
        
        fluidRow(
          
          box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
              
              verticalLayout(
                renderPlot({ Enrichment.plot(object = local.rea.values$dataset.SE, listNames = listname, 
                                             from =       "DiffExpEnrichAnal",
                                             Over_Under = input[[paste0(listname, "-Over_Under")]],
                                             top =        input[[paste0(listname, "-top.over")]],
                                             domain =     input[[paste0(listname, "-domain")]]) }),
                fluidRow(
                  column(3,
                         numericInput(inputId = session$ns(paste0(listname, "-top.over")), label="Top genes :" , value=50 , 
                                      min = 20, max=length(unique(data$Term)), step = 20)),
                  column(4,
                         radioButtons(inputId = session$ns(paste0(listname, "-domain")), label="Domain" ,
                                      choices = unique(data$Domain), selected = unique(data$Domain)[1], inline = FALSE, width = 1.5)),
                  column(5,
                         radioButtons(inputId = session$ns(paste0(listname, "-Over_Under")), label="" ,
                                      choices = unique(data$Decision), selected = unique(data$Decision)[1], inline = FALSE, width = 1.5))
                )),
              hr(),
              verticalLayout(
                
                tags$br(),
                ## DEG result table ###
                DT::renderDataTable({
                  
                  DT::datatable(dplyr::filter(data, Decision == input[[paste0(listname, "-Over_Under")]], 
                                              Domain   == input[[paste0(listname, "-domain")]]),
                                options = list(rownames = TRUE, pageLength = 5, lengthMenu = c(5,10, 15, 20), scrollX = T))
                })
              )
          )
        )
      }
    })
  })
  
  output$AnnotCoExpResults <- renderUI({
    
    if(rea.values[[dataset]]$coExpAnnot == FALSE) return()
    
    #dataset.SE <- FlomicsMultiAssay@ExperimentList[[paste0(dataset,".filtred")]]
    
    
    # check if result is empty
    if(is.null(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]])){
      
      showModal(modalDialog( title = "Error message", "There is no correspondence between the feature lists and the annotation file"))
    }
    validate({ 
      need(!is.null(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]]), message="There is no correspondence between the feature lists and the annotation file") 
    })
    
    # foreach gene list selected (contrast)
    lapply(names(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]]), function(listname) {
      
      # if result is empty
      if(is.null(local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]][[listname]])){
        
        box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
            "There is no correspondence between the feature list and the annotation file !"
        )
        
      }else{
        
        data <- local.rea.values$dataset.SE@metadata$CoExpEnrichAnal[["results"]][[listname]][["Over_Under_Results"]]
        
        fluidRow(
          
          box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = listname, 
              
              verticalLayout(
                renderPlot({ Enrichment.plot(dataset.SE, listNames=listname, from="CoExpEnrichAnal",
                                             Over_Under = input[[paste0(listname, "-Over_Under")]],
                                             top =        input[[paste0(listname, "-top.over")]],
                                             domain =     input[[paste0(listname, "-domain")]]) }),
                fluidRow(
                  column(3,
                         numericInput(inputId = session$ns(paste0(listname, "-top.over")), label="Top genes :" , value=50 , 
                                      min = 20, max=length(unique(data$Term)), step = 20)),
                  column(4,
                         radioButtons(inputId = session$ns(paste0(listname, "-domain")), label="Domain" ,
                                      choices = unique(data$Domain), selected = unique(data$Domain)[1], inline = FALSE, width = 1.5)),
                  column(5,
                         radioButtons(inputId = session$ns(paste0(listname, "-Over_Under")), label="" ,
                                      choices = unique(data$Decision), selected = unique(data$Decision)[1], inline = FALSE, width = 1.5))
                )),
              hr(),
              verticalLayout(
                
                tags$br(),
                ## DEG result table ###
                DT::renderDataTable({
                  
                  DT::datatable(dplyr::filter(data, Decision == input[[paste0(listname, "-Over_Under")]], 
                                              Domain   == input[[paste0(listname, "-domain")]]),
                                options = list(rownames = TRUE, pageLength = 5, lengthMenu = c(5,10, 15, 20), scrollX = T))
                })
              )
          ) 
        )
      }
    })
  })
  
}

