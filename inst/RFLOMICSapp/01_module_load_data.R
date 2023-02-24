####
# LOAD data module
####

LoadOmicsDataUI <- function(id){
  
  
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      # load exp design
      # display tab
      box(width = 12, title = "Load experimental design", status = "warning", solidHeader = TRUE,
          # ExpDesign
          fluidRow(
            column(width = 12,
                   # 1- set project name
                   column(width = 4, textInput(inputId = ns("projectName"), label = "Project name")),
                   # 2- matrix count/abundance input
                   column(width = 8, fileInput(inputId = ns("Experimental.Design.file"), label = "Experimental Design (tsv)",
                                               accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
                   # #metadata
                   # column(width = 4, fileInput(inputId = ns("metadata.file"), label = "metadata (tsv)",
                   #                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")))
            )
          ),
          fluidRow(
            column(width = 12,
                   uiOutput(ns("ExpDesignTable"))
            ),
            column(width = 12,
                   uiOutput(ns("selectSamplesUI"))),
          )
      )
    ),
    fluidRow( uiOutput(outputId = ns("LoadDataUI"))),
    fluidRow( column(width = 6, actionButton(inputId = ns("loadData"),"load Data"))),
    br(),
    
    fluidRow(
      uiOutput(ns("CompletenessUI")),
      uiOutput(ns("summaryMAE"))
    )
    
  )
  
}


LoadOmicsData <- function(input, output, session, rea.values){
  
  local.rea.values <- reactiveValues(plots = FALSE, ExpDesign = NULL, FactorList = NULL)
  
  observe({
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    
    rea.values$contrastMat  <- NULL
    #FlomicsMultiAssay      <<- NULL
    session$userData$FlomicsMultiAssay <- NULL
    
    rea.values$validate.status <- 0
    
  })
  
  
  ##########################################
  # (1) load experimental design
  ##########################################
  
  # as soon as the "design file" has been loaded
  # => check file format
  # => display content
  # => set ref for each factor
  # => set type of factors (bio, batch, meta)
  
  observeEvent(input$Experimental.Design.file, {
    
    # reset
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    local.rea.values$plots <- FALSE
    
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    
    session$userData$FlomicsMultiAssay      <- NULL
    #session$userData$Design <- NULL
    
    print("# 1- Load data ...")
    
    print("#    => Load experimental design...")
    
    # read and check design file
    design.tt <- tryCatch(expr = read_exp_design(file = input$Experimental.Design.file$datapath),
                          error=function(e) e, warning=function(w) w)
    
    if(!is.null(design.tt$message)){
      
      showModal(modalDialog( title = "Error message", design.tt$message)) 
    }
    validate({ need(expr = is.null(design.tt$message), message=design.tt$message) })
    
    ExpDesign.tbl <- design.tt
    
    local.rea.values$ExpDesign <- ExpDesign.tbl
    
    ####### Display tab of exp design  ########
    
    # display table
    output$ExpDesignTable <- renderUI({
      box(width = 12, background = "light-blue", solidHeader = TRUE, collapsible = TRUE, 
          collapsed = TRUE, title = "Overview of experimental design table", 
          
          # DT::renderDataTable( DT::datatable(data = ExpDesign.tbl.affich, filter = 'top',
          #                                    options = list( pageLength = 5, autoWidth = TRUE, dom = 'tp' )))
          #                                    #caption = 'Table 1: This is a simple caption for the table.'
          DT::renderDataTable( DT::datatable(data = ExpDesign.tbl,
                                             options = list( pageLength = 5, autoWidth = TRUE, dom = 'tp' )))
      )
    })
    
    #### order and select modality for each factor
    output$selectSamplesUI <- renderUI({
      
      # condition list
      box(width = 12, background = "green", # Valid colors are: red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
          tags$b("The list and order of conditions :"),
          fluidRow(
            lapply(names(ExpDesign.tbl), function(i) {
              
              column(width = round(12/length(names(ExpDesign.tbl))),
                     # sample list :
                     # orderInput(inputId = session$ns(paste0("selectFactors.", i)),
                     #            label   = tags$span(style="color: black;",i) ,
                     #            items   = levels(as.factor(ExpDesign.tbl[[i]])))
                     selectizeInput(inputId = session$ns(paste0("selectFactors.", i)),
                                    label   = tags$span(style="color: black;",i) ,
                                    choices = levels(as.factor(ExpDesign.tbl[[i]])),
                                    selected= levels(as.factor(ExpDesign.tbl[[i]])), multiple = TRUE,
                                    options = NULL))
            }),
            uiOutput(session$ns("GetdFactorRef")))
      )
    })
    
    ####### set ref and type of each factor ########
    output$GetdFactorRef <- renderUI({
      
      # filter samples
      #ExpDesign.tbl.flt <<- ExpDesign.tbl[input$selectSamples,]
      # filtering per conditions
      for(i in names(ExpDesign.tbl)) {
        
        # if select 1 modality or 0 for a Factor we exclude this factor
        if(length(input[[paste0("selectFactors.", i)]]) > 0){
          ExpDesign.tbl <- dplyr::filter(ExpDesign.tbl, get(i) %in% input[[paste0("selectFactors.", i)]])
          ExpDesign.tbl[[i]] <- factor(ExpDesign.tbl[[i]], levels = input[[paste0("selectFactors.", i)]])
        }
        
        if(length(input[[paste0("selectFactors.", i)]]) <= 1){
          ExpDesign.tbl <- dplyr::select(ExpDesign.tbl, -all_of(i))
        }
      }
      
      local.rea.values$ExpDesign <- ExpDesign.tbl
      
      column(width = 12,
             
             # Construct the form to set the reference factor level
             tags$b("Set the type and the reference fo each design factor :"),
             fluidRow(
               lapply(names(ExpDesign.tbl), function(i){
                 
                 column(width= round(12/length(names(ExpDesign.tbl))),
                        selectInput(inputId = session$ns(paste0("dF.RefLevel.", i)),
                                    label   = tags$span(style="color: black;",i) ,
                                    choices = unique(ExpDesign.tbl[[i]])),
                        
                        radioButtons(session$ns(paste0("dF.Type.", i)), label=NULL , choices = c("Bio","batch", "Meta"),
                                     selected = "Bio", inline = FALSE, width = 2, choiceNames = NULL, choiceValues = NULL))
               })
             )
      )
    })
  })
  
  
  ##########################################
  # (2) load data
  ##########################################
  
  # display interface for load data
  output$LoadDataUI <- renderUI({
    
    box(width = 12, title = "Load omic data", status = "warning",  height = NULL, solidHeader = TRUE,
        #h5("[warning] dataset with omics (type == none) was igored..."),
        fluidRow(
          column(2,
                 # omic type
                 selectInput(inputId = session$ns('omicType1'), label='Omic type', choices = c("None"="none", "RNAseq"="RNAseq", "Proteomics"="proteomics", "Metabolomics"="metabolomics"), selected = "none")
          ),
          column(6,
                 # matrix count/abundance input
                 fileInput(inputId = session$ns("data1"), "Dataset matrix (tsv)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          ),
          # column(4,
          #        # metadata/QC bioinfo
          #        fileInput(inputId = session$ns("metadataQC1"), "QC or metadata (Ex.)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          # ),
          column(3,
                 # dataset Name
                 textInput(inputId = session$ns("DataName1"), label="Dataset name", value="set1")
          )
        ),
        
        uiOutput(   outputId = session$ns("toAddData2")),
        actionButton(inputId = session$ns("addData"),"Add data")
    )
  })
  
  # as soon as the "add data" buttom has been clicked
  # => a new select/file Input was display
  # =>
  addDataNum <- 1
  dataName.vec  <- c()
  observeEvent(input$addData, {
    
    # add input select for new data
    addDataNum <<- addDataNum + 1
    output[[paste("toAddData", addDataNum, sep="")]] <- renderUI({
      list(
        fluidRow(
          column(2,
                 # omic type
                 selectInput(inputId=session$ns(paste0('omicType', addDataNum)), label='Omic type', choices = c("None"="none", "RNAseq"="RNAseq", "Proteomics"="proteomics", "Metabolomics"="metabolomics"), selected = "none")),
          column(6,
                 # matrix count/abundance input
                 fileInput(inputId=session$ns(paste0("data", addDataNum)), "Dataset matrix (tsv)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
          # column(4,
          #        # metadata/QC bioinfo
          #        fileInput(inputId=session$ns(paste0("metadataQC", addDataNum)), "QC or metadata (tsv)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
          column(3,
                 # dataset Name
                 textInput(inputId=session$ns(paste0("DataName", addDataNum)), label="Dataset name", value=paste0("set", as.character(addDataNum))))
        ),
        
        uiOutput(session$ns(paste("toAddData",addDataNum + 1,sep="")))
      )
    })
  })
  
  # as soon as the "load" buttom has been clicked
  # => create ExpDesign object
  # => create flomicsMultiAssay object
  # => upsetR
  observeEvent(input$loadData, {
    
    #print(paste0("loadData ", input$loadData))
    ### load Design
    
    # reset objects and UI
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    rea.values$Contrasts.Sel <- NULL
    
    session$userData$FlomicsMultiAssay        <- NULL
    #session$userData$Design   <- NULL
    
    local.rea.values$plots <- FALSE
    # #updateTabItems(session, "tabs", selected = "importData")
    
    ### check project name
    if(input$projectName == ""){
      showModal(modalDialog(title = "Error message", "project name is required"))
    }
    validate({ need(input$projectName != "", message="project name is required") })
    
    
    ### check Experimental Design
    if(is.null(input$Experimental.Design.file)){
      showModal(modalDialog(title = "Error message", "Experimental Design is required"))
    }
    validate({ need(! is.null(input$Experimental.Design.file), message="Set a name") })
    
    
    ExpDesign.tbl <- local.rea.values$ExpDesign
    
    ### Get the Type and ref of the factors that the users enter in the form
    dF.Type.dFac <- vector()
    dF.List.Name <- vector()
    dF.List.ref  <- vector()
    
    rea.values$validate.status <- 0
    
    #factor.names <- stringr::str_remove(stringr::str_subset(string = names(input), pattern = "dF.Type."), pattern = "dF.Type.")
    
    for(dFac in names(ExpDesign.tbl)){
      
      # list of type of factors (bio or batch)
      dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
      # list of level reference of factors
      dF.List.ref[dFac]  <- input[[paste0("dF.RefLevel.",dFac)]]
    }
    
    # check number of factor bio
    if(! length(stringr::str_subset(dF.Type.dFac, "Bio")) %in% 1:3){
      showModal(modalDialog(title = "Error message", "1 to 3 bio factor(s)")) }
    
    # check number of factor batch
    if(! length(stringr::str_subset(dF.Type.dFac, "batch")) %in% 1:2){
      showModal(modalDialog(title = "Error message", "at least 1 batch factor (max = 2)"))}
    
    validate({ need((length(stringr::str_subset(dF.Type.dFac, "Bio"  )) %in% 1:3) &
                      (length(stringr::str_subset(dF.Type.dFac, "batch")) %in% 1:2), message="") })
    
    #### load omics data
    ####
    print("#    => Load omics data...")
    inputs <- list()
    #### list of omic data laoded from interface
    rea.values$validate.status <- 0
    dataName.vec <- c()
    for (k in 1:addDataNum){
      
      if(input[[paste0("omicType", k)]] != "none"){
        
        ### omics type ###
        omicType <- input[[paste0("omicType", k)]]
        
        ### dataset name ### 
        # => check presence of dataname
        dataName.tmp <- gsub("[[:space:]]", "", input[[paste0("DataName", k)]])
        
        if(dataName.tmp == ""){
          showModal(modalDialog( title = "Error message", "Dataset names is required : dataset ", k ))
          rea.values$validate.status <- 1
        }
        dataName <- paste0(omicType, ".", dataName.tmp)
        dataName.vec <- c(dataName.vec, dataName)
        
        # => check duplicat dataset name
        if(any(duplicated(dataName.vec)) == TRUE){
          showModal(modalDialog( title = "Error message", "Dataset names must be unique : dataset ", (1:addDataNum)[duplicated(dataName.vec)] ))
          rea.values$validate.status <- 1
        }
        
        #### omics dataset
        # => check omics data
        if(is.null(input[[paste0("data", k)]])){
          showModal(modalDialog( title = "Error message", "omics dataset is required : dataset ", k ))
          rea.values$validate.status <- 1
        }
        validate({ need(expr = !is.null(input[[paste0("data", k)]]), message="error") })
        
        # => read data matrix
        dataFile <- input[[paste0("data", k)]]
        data.mat.tt <- tryCatch(read_omics_data(file = dataFile$datapath), error=function(e) e, warning=function(w) w)
        
        if(!is.null(data.mat.tt$message)){
          
          showModal(modalDialog( title = "Error message", data.mat.tt$message)) 
          rea.values$validate.status <- 1
        }
        validate({ need(expr = is.null(data.mat.tt$message), message=data.mat.tt$message) })
        
        data.mat <- data.mat.tt
        inputs[[dataName]][["omicType"]] <- omicType
        inputs[[dataName]][["data"]] <- data.mat
        
        #### read meta data file
        # if(!is.null(input[[paste0("metadataQC", k)]])){
        #   print("# ... metadata QC...")
        #   qcFile <- input[[paste0("metadataQC", k)]]
        #   QCmat  <- read_metaData(qcFile$datapath)
        #   
        #   validate({ need(!is.null(QCmat), message="") })
        #   
        #   inputs[[dataName]][["meta"]] <- QCmat
        #   #inputs[[dataName]][["qcFile"]] <- input[[paste0("metadataQC", k)]]$datapath ####
        # }
        inputs[[dataName]][["meta"]]     <- NULL
        
        validate({ need(rea.values$validate.status == 0, message="error") })
      }
    }
    
    # FlomicsMultiAssay.try <- RFLOMICS::try_rflomics(FlomicsMultiAssay.constructor(inputs = inputs,
    #                                                                     ExpDesign = ExpDesign.tbl,
    #                                                                     refList = dF.List.ref,
    #                                                                     typeList = dF.Type.dFac,
    #                                                                     projectName=input$projectName))
    # 
    # if(!is.null(FlomicsMultiAssay.try[["error"]])) {
    #   showModal(modalDialog( title = "Error message", as.character(FlomicsMultiAssay.try[["error"]])))
    # }
    # validate({ need(is.null(FlomicsMultiAssay.try[["error"]]), message="error") })
    
    
    FlomicsMultiAssay.try <- tryCatch(FlomicsMultiAssay.constructor(inputs = inputs,
                                                                    ExpDesign = ExpDesign.tbl,
                                                                    refList = dF.List.ref,
                                                                    typeList = dF.Type.dFac,
                                                                    projectName=input$projectName),
                                      error=function(e) e, warning=function(w) w)
    
    if(!is.null(FlomicsMultiAssay.try$message)) {
      showModal(modalDialog( title = "Error message", FlomicsMultiAssay.try$message))
      rea.values$validate.status <- 1
    }
    validate({ need(is.null(FlomicsMultiAssay.try$message), message="error") })
    
    session$userData$FlomicsMultiAssay <- FlomicsMultiAssay.try
    
    # FlomicsMultiAssay[, complete.cases(FlomicsMultiAssay), ]
    # colData(FlomicsMultiAssay)[!complete.cases(FlomicsMultiAssay),]
    # intersectColumns(FlomicsMultiAssay)
    
    # => completeness
    #### check experimental design : experimental design must be a complete and balanced.
    
    print(paste0("#    => Design Completeness Check..."))
    
    local.rea.values$completeCheckRes <- CheckExpDesignCompleteness(object = session$userData$FlomicsMultiAssay)
    
    # plot : ok
    local.rea.values$plots <- TRUE
    
    if(!is.null(local.rea.values$completeCheckRes[["error"]])){
      #rea.values$loadData <- FALSE
      #rea.values$model    <- FALSE
      showModal(modalDialog(title = "Error message", local.rea.values$completeCheckRes[["error"]]))
      rea.values$validate.status <- 1
    }
    
    # continue only if message is true or warning
    validate({ need(is.null(local.rea.values$completeCheckRes[["error"]]) ,message="") })
    
    # 
    rea.values$datasetList <- session$userData$FlomicsMultiAssay@metadata$omicList
    
    rea.values$loadData <- TRUE
    #rea.values$model    <- TRUE
    
  }, ignoreInit = TRUE)
  
  ##########################################
  # (3) check data
  ##########################################
  
  # completeness check
  output$CompletenessUI <- renderUI({
    
    if (local.rea.values$plots == FALSE) return()
    #if (rea.values$loadData == FALSE) return()
    
    print(paste0("#    => Completeness plot..."))
    
    box( width = 6,  status = "warning", title = "Completeness check", solidHeader = TRUE,
         
         # plot of count per condition
         renderPlot(
           isolate({ local.rea.values$completeCheckRes[["plot"]] })
         ),
         hr(),
         tags$i("You **must** have a **complete design** (i.e. all possible combinations of factor's level).
                     **Balanced design** (presence of the same number of replicats for all
                     possible combinations) is not required  but advised.
                     You **must** also have at least one biological factor and 2 replicats")
    )
  })
  
  # upset of all data
  output$summaryMAE <- renderUI({
    
    if (local.rea.values$plots == FALSE) return()
    #if (rea.values$loadData == FALSE) return()
    
    print(paste0("#    => overview plot..."))
    
    box(width = 6, status = "warning", title = "Dataset(s) overview", solidHeader = TRUE,
        renderPlot( isolate({ RFLOMICS::Datasets_overview_plot(session$userData$FlomicsMultiAssay) })),
    )
    #}
  })
  
  return(input)
}


######### FUNCTION ##########

read_metaData <- function(input.file){
  
  res <- RFLOMICS::try_rflomics(read.table(input.file, header = TRUE,row.names = 1, sep = "\t"))
  
  # If an error occured
  if(!is.null(res[["error"]]))
  {
    showModal(modalDialog( title = "Error message", as.character(res[["error"]])))
    return(NULL)
  }
  else{
    
    if(!is.null(res$primary) | !is.null(res$colname)){
      
      showModal(modalDialog( title = "Error message", "metadata file must contient at least 2 colomns named : primary (?) and colname(?)"))
      return(NULL)
      
    }
  }
  
  return(res$value)
}
