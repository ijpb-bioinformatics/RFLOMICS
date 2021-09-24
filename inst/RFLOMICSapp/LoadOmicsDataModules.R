


LoadOmicsDataUI <- function(id){


  ns <- NS(id)

  tagList(
    
    # set project name
    fluidRow(
        box(width = 12, status = "warning", solidHeader = TRUE, 
            title = "Load omics data and experimental design",
            column(4, 
                   textInput(inputId = ns("projectName"), label = "Project name")))),
    
    
    fluidRow(
        # load exp design
        # display tab
        box(width = 6, title = "Load experimental design", status = "warning",
            
            # matrix count/abundance input
            fileInput(inputId = ns("Experimental.Design.file"), label = "Import matrix of Experimental Design (txt)",
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
            tags$br(),
            # table visualisation
            dataTableOutput(ns("ExpDesignTable")),
            tags$br(),
            #uiOutput(ns("ExpDesignTable")),
            uiOutput(ns("GetdFactorRef")),
            tags$br(),
            tags$br(),
            uiOutput(ns("CompletenessUI")),
            ),
        
        # load exp design
        # display tab
        uiOutput(ns("LoadDataUI")),
        uiOutput(ns("UpsetMAE"))
    )
  )

}


LoadOmicsData <- function(input, output, session){

    # as soon as the "design file" has been loaded
    # => check file format
    # => display content
    # => set ref for each factor
    # => set type of factors (bio, batch, meta)
    observeEvent(input$Experimental.Design.file, {
    
        # read and check design file
        ExpDesign.tbl <- read_ExpDesign(input.file = input$Experimental.Design.file$datapath)
    
        validate({ need(!is.null(ExpDesign.tbl), message="") })
        
        print("# 1- Load experimental design...")
        
        ####### Display tab of exp design  ########
        if(dim(ExpDesign.tbl)[2] > 5){ ExpDesign.tbl.affich <- ExpDesign.tbl[,1:5] }else{ ExpDesign.tbl.affich <- ExpDesign.tbl }
        
        # output$ExpDesignTable <- renderUI({
        #   box(12, DT::renderDataTable( DT::datatable(ExpDesign.tbl, options = list(pageLength = 5, searching = FALSE)))   ) })
        
        output$ExpDesignTable <- renderDataTable({ DT::datatable(ExpDesign.tbl.affich, options = list(pageLength = 5, searching = FALSE) ) })
    
        ####### Set up Design model ########
        output$GetdFactorRef <- renderUI({
          
          tagList(
            # Construct the form to set the reference factor level
            h4("Select the level of reference fo each design factor"),
            fluidRow(
              lapply(names(ExpDesign.tbl.affich), function(i) {
                box(width=3,
                    selectInput(session$ns(paste0("dF.RefLevel.", i)), i, choices = levels(as.factor(ExpDesign.tbl.affich[[i]])), selectize=FALSE, size=5))})),
            
            # Construct the form to set the type of the factor (either biological or batch)
            h4("Select the type of the design factor"),
            fluidRow(
              lapply(names(ExpDesign.tbl.affich), function(i) {
                box(width=3,
                    radioButtons(session$ns(paste0("dF.Type.", i)), label=NULL , choices = c("Bio","batch", "Meta"), selected = "Bio", inline = FALSE, width = 2, choiceNames = NULL, choiceValues = NULL))})),
            actionButton(session$ns("ValidF"),"Valid factor set up")
          )
        })
  
    })
 
    ##########################################
    # Part2 : Define the Experimental design:
    #         -> load experimental plan
    #         -> the level of ref for each factor
    ##########################################
     
    # as soon as the "valid factor" buttom has been clicked
    # => 
    # =>
    observeEvent(input$ValidF, {
      
      ### check project name
      if(input$projectName == ""){
        showModal(modalDialog(title = "Error message", "project name is required"))
      }
      validate({ need(input$projectName != "", message="project name is required") })
      
      
      ### Experimental Design
      if(is.null(input$Experimental.Design.file)){
        showModal(modalDialog(title = "Error message", "Experimental Design is required"))
      }
      validate({ need(! is.null(input$Experimental.Design.file), message="Set a name") })
      
      ExpDesign.tbl <- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1, sep = "\t")
      
      ### Get the Type and ref of the factors that the users enter in the form
      dF.Type.dFac<-vector()
      dF.List.Name<-vector()
      dF.List.ref <-vector()
      
      validate.status <<- 0
      
      #factor.names <- stringr::str_remove(stringr::str_subset(string = names(input), pattern = "dF.Type."), pattern = "dF.Type.")      

      for(dFac in names(ExpDesign.tbl)){
        
          # list of type of factors (bio or batch)
          dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
          # list of level reference of factors
          dF.List.ref[dFac]  <- input[[paste0("dF.RefLevel.",dFac)]]
          }
      
      if(! length(stringr::str_subset(dF.Type.dFac, "Bio")) %in% 1:3){
        showModal(modalDialog(title = "Error message", "1 to 3 bio factor(s)")) }
      
      if(length(stringr::str_subset(dF.Type.dFac, "batch")) == 0){
        showModal(modalDialog(title = "Error message", "at least 1 batch factor"))}     
      
      validate({ need((length(stringr::str_subset(dF.Type.dFac, "Bio"  )) %in% 1:3) & 
                      (length(stringr::str_subset(dF.Type.dFac, "batch")) != 0), message="") })
      
      ## create Design object
      Design <<- ExpDesign.constructor(ExpDesign = ExpDesign.tbl, refList = dF.List.ref, typeList = dF.Type.dFac)
      
      
      # => completeness
      #### check experimental design : experimental design must be a complete and balanced. 
      
      print(paste0("# 3- Design Completeness Check..."))
      
      completeCheckRes <- CheckExpDesignCompleteness(Design)
      if(completeCheckRes[["message"]][1] == "false"){
        showModal(modalDialog(title = "Error message", completeCheckRes[["message"]][2]))
        
      }
      
      output$CompletenessUI <- renderUI({ 
        
        box( width = 12, title = "Completeness",  status = "warning",
             
             #print message
             renderText( completeCheckRes[["message"]][2] ),
             hr(),
             # plot of count per condition
             renderPlot( plotExperimentalDesign(completeCheckRes[["count"]] )),
             hr(),
             tags$i("You **must** have a **complete design** (i.e. all possible combinations of factor's level).
                     **Balanced design** (presence of the same number of replicats for all
                     possible combinations) is not required  but advised.
                     You **must** also have at least one biological factor and 2 replicats")
        )
      })

      # continue only if message is true or warning
      validate({ need(completeCheckRes[["message"]][1] != "false" ,message="ok") })
      
      
      ### display interface for load data
      output$LoadDataUI <- renderUI({
      
        box(width = 6, title = "Load omic data", status = "warning",  height = NULL,
            #h5("[warning] dataset with omics (type == none) was igored..."),
            fluidRow(
              column(2,
                     # omic type
                     selectInput(inputId = session$ns('omicType1'), label='Omics', choices = c("None"="none", "RNAseq"="RNAseq", "Proteomics"="proteomics", "Metabolomics"="metabolomics"), selected = "none")
              ),
              column(4,
                     # matrix count/abundance input
                     fileInput(inputId = session$ns("data1"), "omics count/abundance (Ex.)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
              ),
              column(4,
                     # metadata/QC bioinfo
                     fileInput(inputId = session$ns("metadataQC1"), "QC or metadata (Ex.)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
              ),
              column(2,
                     # dataset Name
                     textInput(inputId = session$ns("DataName1"), label="Dataset name", value="set1")
              )
            ),
            # fluidRow(
            #          uiOutput(session$ns("CompletenessUI1"))
            # ),
            uiOutput(   outputId = session$ns("toAddData2")),
            actionButton(inputId = session$ns("addData"),"Add data"),
            actionButton(inputId = session$ns("loadData"),"load")
        )
      })
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
                     selectInput(inputId=session$ns(paste0('omicType', addDataNum)), label='Omics', choices = c("None"="none", "RNAseq"="RNAseq", "Proteomics"="proteomics", "Metabolomics"="metabolomics"), selected = "none")),
              column(4,
                     # matrix count/abundance input
                     fileInput(inputId=session$ns(paste0("data", addDataNum)), "Feature count/abundance (txt)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
              column(4,
                     # metadata/QC bioinfo
                     fileInput(inputId=session$ns(paste0("metadataQC", addDataNum)), "QC or metadata (txt)", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
              column(2,
                     # dataset Name
                     textInput(inputId=session$ns(paste0("DataName", addDataNum)), label="Dataset name", value=paste0("set", as.character(addDataNum))))
            ),
            # fluidRow(
            #   column(12,
            #          uiOutput(outputId = session$ns(paste0("CompletenessUI", addDataNum))))
            #   ),
            
            uiOutput(session$ns(paste("toAddData",addDataNum + 1,sep="")))
          )
        })
    })
    
    
    # as soon as the "load" buttom has been clicked
    # => create ExpDesign object
    # => create flomicsMultiAssay object
    # => upsetR
    observeEvent(input$loadData, {

      #### load data
      ####
      print("# 2- Load omics data...")
      inputs <- list()
      omicList <- list()
      #### list of omic data laoded from interface
      validate.status <<- 0
      dataName.vec <- c()
      for (k in 1:addDataNum){

        if(input[[paste0("omicType", k)]] != "none"){

          omicType <- input[[paste0("omicType", k)]]

          dataName <- paste0(omicType, ".", gsub("[[:space:]]", "", input[[paste0("DataName", k)]]))

          
          
          # check presence of dataname
          if(dataName == ""){
            showModal(modalDialog( title = "Error message", "Dataset names is required : dataset ", k ))
            validate.status <<- 1
          }

          # check duplicat dataset name
          if(any(duplicated(omicList[[omicType]])) == TRUE){
            showModal(modalDialog( title = "Error message", "Dataset names must be unique : dataset ", (1:addDataNum)[duplicated(dataName.vec)] ))
            validate.status <<- 1
          }

          #### read omic data file
          # => check omics data
          if(is.null(input[[paste0("data", k)]])){
            showModal(modalDialog( title = "Error message", "omics counts/abundances matrix is required : dataset ", k ))
            validate.status <<- 1
          }
          dataFile <- input[[paste0("data", k)]]
          data.mat <- read_omicsData(dataFile$datapath)
          
          validate({ need(!is.null(data.mat), message="") })
          
          inputs[[dataName]][["data"]]     <- data.mat
          inputs[[dataName]][["omicType"]] <- omicType
          inputs[[dataName]][["meta"]]     <- NULL
          #inputs[[dataName]][["dataFile"]] <- input[[paste0("data", k)]]$datapath 
          
          #### read meta data file
          if(!is.null(input[[paste0("metadataQC", k)]])){
            print("# ... metadata QC...")
            qcFile <- input[[paste0("metadataQC", k)]]
            QCmat  <- read_metaData(qcFile$datapath)
            
            validate({ need(!is.null(QCmat), message="") })
            
            inputs[[dataName]][["meta"]] <- QCmat
            #inputs[[dataName]][["qcFile"]] <- input[[paste0("metadataQC", k)]]$datapath ####
          }

          validate({ need(validate.status == 0, message="error") })
        }
      }

      FlomicsMultiAssay.try <- try_rflomics(FlomicsMultiAssay.constructor(inputs = inputs, Design=Design, projectName=input$projectName))

      if(!is.null(FlomicsMultiAssay.try[["error"]])) {
        showModal(modalDialog( title = "Error message", as.character(FlomicsMultiAssay.try[["error"]])))
      }
      validate({ need(is.null(FlomicsMultiAssay.try[["error"]]), message="error") })
      
      FlomicsMultiAssay <<- FlomicsMultiAssay.try[["value"]]
      
      # FlomicsMultiAssay[, complete.cases(FlomicsMultiAssay), ]
      # colData(FlomicsMultiAssay)[!complete.cases(FlomicsMultiAssay),]
      # intersectColumns(FlomicsMultiAssay)   
      

        
      
      # upset of all data
      output$UpsetMAE <- renderUI({
        
        if (length(names(FlomicsMultiAssay)) >1){
          box(width = 6, status = "warning", renderPlot(upsetSamples(FlomicsMultiAssay)))}
      })

    })

  return(input)

}


######### FUNCTION ##########

read_ExpDesign <- function(input.file){
  
  res <- try_rflomics(read.table(input.file, header = TRUE,row.names = 1, sep = "\t"))
  
  # If an error occured
  if(!is.null(res[["error"]]))
  {
    showModal(modalDialog( title = "Error message", as.character(res[["error"]])))
    return(NULL)
  }
  else{
    
    if(length(dim(res$value)[2]) > 5){ value <- res$value[,1:5] }else{ value <- res$value }
    
    # check nbr of modality of the 5th fist columns 
    index <- sapply(names(value), function(x){ if(length(unique(value[[x]]))>10){ FALSE }else{ TRUE } })
    F.mod <- names(value)[index]
    
    ratio <- length(F.mod)/length(names(value))
    
    if(ratio != 1)
    {
      showModal(modalDialog( title = "Warning message", "The select input contains a large number of options"))
    }
    
    return(res$value)
  }
}


read_omicsData <- function(input.file){
  
  res <- try_rflomics(read.table(input.file, header = TRUE,row.names = 1, sep = "\t"))
  
  # If an error occured
  if(!is.null(res[["error"]]))
  {
    showModal(modalDialog( title = "Error message", as.character(res[["error"]])))
    return(NULL)
  }
    
  return(res$value)
}

read_metaData <- function(input.file){
  
  res <- try_rflomics(read.table(input.file, header = TRUE,row.names = 1, sep = "\t"))
  
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
