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
                   column(width = 3, textInput(inputId = ns("projectName"), label = "Project name")),
                   # 2- matrix count/abundance input
                   column(width = 7, fileInput(inputId = ns("Experimental.Design.file"), label = "Experimental Design (tsv)",
                                               accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv"))),
                   column(width = 2, actionButton(inputId = ns("loadEx_metadata"), label = "Load Example Metadata", icon = shiny::icon("file-import"))),
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
    fluidRow( column(width = 2, actionButton(inputId = ns("loadData"),"load Data")),
              column(width = 2, actionButton(inputId = ns("loadExData"),"load Example Data", icon = shiny::icon("file-import")))),
    br(),
    
    fluidRow(
      uiOutput(ns("CompletenessUI")),
      uiOutput(ns("summaryMAE"))
    )
    
  )
  
}

LoadOmicsData <- function(input, output, session, rea.values){
  
  local.rea.values <- reactiveValues(plots = FALSE, ExpDesign = NULL, FactorList = NULL, dataPath = NULL, listInputs = NULL)
  
  observe({
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    
    rea.values$contrastMat  <- NULL
    session$userData$FlomicsMultiAssay <- NULL
    
    rea.values$validate.status <- 0
    
  })
  
  
  # ---- Datapath depending on the button ----
  # load user own metadata file
  observeEvent(input$Experimental.Design.file, {
    local.rea.values$dataPath <- NULL
    local.rea.values$dataPath <- input$Experimental.Design.file$datapath
  })
  
  # load example metadata file
  observeEvent(input$loadEx_metadata, {
    local.rea.values$dataPath <- NULL
    local.rea.values$dataPath <- paste0(system.file(package = "RFLOMICS"), 
                                        "/ExamplesFiles/ecoseed/condition.txt")
    
  })
  
  # ---- Load experimental design ----
  ##########################################
  # (1) load experimental design
  ##########################################
  
  # as soon as the "design file" has been loaded
  # => check file format
  # => display content
  # => set ref for each factor
  # => set type of factors (bio, batch, meta)
  
  observeEvent(local.rea.values$dataPath, {
    
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
    
    if (local.rea.values$dataPath == paste0(system.file(package = "RFLOMICS"), 
                                            "/ExamplesFiles/ecoseed/condition.txt")) {
      print("#    => Load example experiment design...")
    }else print("#    => Load experimental design...")
    
    # read and check design file
    design.tt <- tryCatch(expr =  read_exp_design(file = local.rea.values$dataPath),
                          error = function(e) e, warning = function(w) w)
    
    
    if(!is.null(design.tt$message)){
      
      showModal(modalDialog( title = "Error message", design.tt$message)) 
    }
    validate({ need(expr = is.null(design.tt$message), message = design.tt$message) })
    
    ExpDesign.tbl <- design.tt
    
    local.rea.values$ExpDesign <- ExpDesign.tbl
    
    ####### Display tab of exp design  ########
    
    # display table
    output$ExpDesignTable <- renderUI({
      box(width = 12, background = "light-blue", solidHeader = TRUE, collapsible = TRUE, 
          collapsed = TRUE, title = "Overview of experimental design table", 
          DT::renderDataTable( DT::datatable(data = ExpDesign.tbl,
                                             options = list( pageLength = 5, autoWidth = TRUE, dom = 'tp' )))
      )
    })
    
    #### order and select modality for each factor
    output$selectSamplesUI <- renderUI({
      
      # condition list
      box(width = 12, background = "green", # Valid colors are: red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
          tags$b("The list and order of conditions:"),
          fluidRow(
            lapply(names(ExpDesign.tbl), function(i) {
              
              column(width = round(12/length(names(ExpDesign.tbl))),
                     selectizeInput(inputId = session$ns(paste0("selectFactors.", i)),
                                    label   = tags$span(style = "color: black;",i) ,
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
      # filtering per conditions
      for(i in names(ExpDesign.tbl)) {
        
        # if select 1 modality or 0 for a Factor we exclude this factor
        if(length(input[[paste0("selectFactors.", i)]]) > 0){
          ExpDesign.tbl <- dplyr::filter(ExpDesign.tbl, get(i) %in% input[[paste0("selectFactors.", i)]])
          ExpDesign.tbl[[i]] <- factor(ExpDesign.tbl[[i]], levels = input[[paste0("selectFactors.", i)]])
        }
        
        if(length(input[[paste0("selectFactors.", i)]]) <= 1){
          ExpDesign.tbl <- dplyr::select(ExpDesign.tbl, -tidyselect::all_of(i))
        }
      }
      
      local.rea.values$ExpDesign <- ExpDesign.tbl
      
      column(width = 12,
             
             # Construct the form to set the reference factor level
             tags$b("Set the type and the reference fo each design factor:"),
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
  
  # ---- Load Data ----
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
          column(3,
                 # dataset Name
                 textInput(inputId = session$ns("DataName1"), label="Dataset name", value="set1")
          )
        ),
        uiOutput(   outputId = session$ns("toAddData2")),
        actionButton(inputId = session$ns("addData"),   "Add data")
    )
  })
  
  
  # ---- observe Event add Data ----
  # as soon as the "add data" button has been clicked
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
                 selectInput(inputId = session$ns(paste0('omicType', addDataNum)), label = 'Omic type', choices = c("None" = "none", "RNAseq" = "RNAseq", "Proteomics" = "proteomics", "Metabolomics" = "metabolomics"), selected = "none")),
          column(6,
                 # matrix count/abundance input
                 fileInput(inputId = session$ns(paste0("data", addDataNum)), "Dataset matrix (tsv)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
          column(3,
                 # dataset Name
                 textInput(inputId = session$ns(paste0("DataName", addDataNum)), label = "Dataset name", value=paste0("set", as.character(addDataNum))))
        ),
        
        uiOutput(session$ns(paste("toAddData", addDataNum + 1, sep = "")))
      )
    })
  })
  
  # ---- observeEvent loadData | loadExData ----
  
  observeEvent(input$loadExData, {
    
    local.rea.values$listInputs <- NULL
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    rea.values$Contrasts.Sel <- NULL
    
    session$userData$FlomicsMultiAssay <- NULL
    
    local.rea.values$plots <- FALSE
    
    print("#    => Load example omics data...")
    
    inputs <- list()
    rea.values$validate.status <- 0
    
    # RNASeq
    data.mat.tt <- tryCatch( read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")), 
                            error = function(e) e, warning = function(w) w)
    dataName <- "RNAseq.set1"
    inputs[[dataName]] <- list("omicType" = "RNAseq", "data" = data.mat.tt, "meta" = NULL) 
    
    # Metabolomic
    data.mat.tt <- tryCatch( read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
                            error = function(e) e, warning = function(w) w)
    dataName <- "metabolomics.set2"
    inputs[[dataName]] <- list("omicType" = "metabolomics", "data" = data.mat.tt, "meta" = NULL) 
    
    # proteomics
    data.mat.tt <- tryCatch( read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")), 
                            error = function(e) e, warning = function(w) w)
    dataName <- "proteomics.set3"
    inputs[[dataName]] <- list("omicType" = "proteomics", "data" = data.mat.tt, "meta" = NULL) 
    
    local.rea.values$listInputs <- inputs
    
    
  })
  
  observeEvent(input$loadData, {
    
    local.rea.values$listInputs <- NULL
    
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    rea.values$Contrasts.Sel <- NULL
    
    session$userData$FlomicsMultiAssay <- NULL
    
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
          showModal(modalDialog( title = "Error message", "Dataset names is required: dataset ", k ))
          rea.values$validate.status <- 1
        }
        dataName <- paste0(omicType, ".", dataName.tmp)
        dataName.vec <- c(dataName.vec, dataName)
        
        # => check duplicat dataset name
        if(any(duplicated(dataName.vec)) == TRUE){
          showModal(modalDialog( title = "Error message", "Dataset names must be unique: dataset ", (1:addDataNum)[duplicated(dataName.vec)] ))
          rea.values$validate.status <- 1
        }
        
        #### omics dataset
        # => check omics data
        if(is.null(input[[paste0("data", k)]])){
          showModal(modalDialog( title = "Error message", "omics dataset is required: dataset ", k ))
          rea.values$validate.status <- 1
        }
        validate({ need(expr = !is.null(input[[paste0("data", k)]]), message="error") })
        
        # => read data matrix
        dataFile <- input[[paste0("data", k)]]
        data.mat.tt <- tryCatch( read_omics_data(file = dataFile$datapath), error=function(e) e, warning=function(w) w)
        
        if(!is.null(data.mat.tt$message)){
          
          showModal(modalDialog( title = "Error message", data.mat.tt$message)) 
          rea.values$validate.status <- 1
        }
        validate({ need(expr = is.null(data.mat.tt$message), message=data.mat.tt$message) })
        
        data.mat <- data.mat.tt
        inputs[[dataName]][["omicType"]] <- omicType
        inputs[[dataName]][["data"]] <- data.mat
        inputs[[dataName]][["meta"]]     <- NULL
        
        validate({ need(rea.values$validate.status == 0, message = "error") })
      }
    }
    
    local.rea.values$listInputs <- inputs
    
  })
  
  # ---- observe Event load Data ----
  # as soon as the "load" buttom has been clicked
  # => create ExpDesign object
  # => create flomicsMultiAssay object
  # => upsetR
  observeEvent(local.rea.values$listInputs, {

    ### load Design
    
    # reset objects and UI
    rea.values$loadData <- FALSE
    rea.values$model    <- FALSE
    rea.values$analysis <- FALSE
    rea.values$Contrasts.Sel <- NULL
    rea.values$datasetList   <- NULL
    rea.values$datasetDiff   <- NULL
    rea.values$Contrasts.Sel <- NULL
    
    session$userData$FlomicsMultiAssay <- NULL
    #session$userData$Design   <- NULL
    
    local.rea.values$plots <- FALSE
    # #updateTabItems(session, "tabs", selected = "importData")
    
    ### check listInputs # no reason to check for null?
    if (is.null(local.rea.values$listInputs)) {
      showModal(modalDialog(title = "Error message", "Please load data"))
    }
    validate({ need(!is.null(local.rea.values$listInputs), message = "Please load data") })
    
    if (length(local.rea.values$listInputs) == 0) {
      showModal(modalDialog(title = "Error message", "Please load data"))
    }
    validate({ need(length(local.rea.values$listInputs) > 0, message = "Please load data") })
  
    ### check project name
    if (input$projectName == "") {
      showModal(modalDialog(title = "Error message", "project name is required"))
    }
    validate({ need(input$projectName != "", message="project name is required") })
    
    if (is.null(local.rea.values$dataPath)) {
      showModal(modalDialog(title = "Error message", "Experimental Design is required"))
    }
    validate({ need(! is.null(local.rea.values$dataPath), message = "Experimental Design is required") })
    
    ExpDesign.tbl <- local.rea.values$ExpDesign
    
    ### Get the Type and ref of the factors that the users enter in the form
    dF.Type.dFac <- vector()
    dF.List.Name <- vector()
    dF.List.ref  <- vector()
    
    rea.values$validate.status <- 0
    
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
    inputs <- list()
    #### list of omic data laoded from interface
    rea.values$validate.status <- 0
   
    FlomicsMultiAssay.try <- tryCatch( FlomicsMultiAssay.constructor(inputs = local.rea.values$listInputs,
                                                                              ExpDesign   = ExpDesign.tbl,
                                                                              refList     = dF.List.ref,
                                                                              typeList    = dF.Type.dFac,
                                                                              projectName = input$projectName),
                                      error = function(e) e, warning = function(w) w)
    
    if(!is.null(FlomicsMultiAssay.try$message)) {
      showModal(modalDialog( title = "Error message", FlomicsMultiAssay.try$message))
      rea.values$validate.status <- 1
    }
    validate({ need(is.null(FlomicsMultiAssay.try$message), message="error") })
    
    session$userData$FlomicsMultiAssay <- FlomicsMultiAssay.try
    
    # => completeness
    #### check experimental design : experimental design must be a complete and balanced.
    
    print(paste0("#    => Design Completeness Check..."))
    
    local.rea.values$completeCheckRes <-  CheckExpDesignCompleteness(object = session$userData$FlomicsMultiAssay)
    
    # plot : ok
    local.rea.values$plots <- TRUE
    
    if (!is.null(local.rea.values$completeCheckRes[["error"]])){
      showModal(modalDialog(title = "Error message", local.rea.values$completeCheckRes[["error"]]))
      rea.values$validate.status <- 1
    }
    
    # continue only if message is true or warning
    validate({ need(is.null(local.rea.values$completeCheckRes[["error"]]), message = local.rea.values$completeCheckRes[["error"]]) })
    
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
    
    print(paste0("#    => Completeness plot..."))
    
    box( width = 6,  status = "warning", title = "Completeness check", solidHeader = TRUE,
         
         # plot of count per condition
         renderPlot(
           isolate({ local.rea.values$completeCheckRes[["plot"]] })
         ),
         hr(),
         tags$div(
           HTML("<em>You <b>must</b> have a <b>complete design</b> (i.e. all possible combinations of factor's levels).
                     <b>Balanced design</b> (presence of the same number of replicates for all
                     possible combinations) is not required  but advised.
                     You <b>must</b> also have at least one biological factor and 2 replicates</em>")
         )
    )
  })
  
  # upset of all data
  output$summaryMAE <- renderUI({
    
    if (local.rea.values$plots == FALSE) return()
    
    print(paste0("#    => overview plot..."))
    
    box(width = 6, status = "warning", title = "Dataset(s) overview", solidHeader = TRUE,
        renderPlot( isolate({  Datasets_overview_plot(session$userData$FlomicsMultiAssay) })),
    )
    #}
  })
  
  return(input)
}


######### FUNCTION ##########

read_metaData <- function(input.file){
  
  res <-  try_rflomics(read.table(input.file, header = TRUE,row.names = 1, sep = "\t"))
  
  # If an error occured
  if(!is.null(res[["error"]]))
  {
    showModal(modalDialog( title = "Error message", as.character(res[["error"]])))
    return(NULL)
  }
  else{
    
    if(!is.null(res$primary) | !is.null(res$colname)){
      
      showModal(modalDialog( title = "Error message", "metadata file must contient at least 2 colomns named: primary (?) and colname(?)"))
      return(NULL)
      
    }
  }
  
  return(res$value)
}
