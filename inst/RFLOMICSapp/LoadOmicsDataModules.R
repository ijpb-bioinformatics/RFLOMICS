


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
                  column(width = 4, fileInput(inputId = ns("Experimental.Design.file"), label = "Experimental Design (txt)",
                                              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
                  #metadata
                  column(width = 4, fileInput(inputId = ns("metadata.file"), label = "metadata (txt)",
                                              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")))
                )
              ),
              fluidRow(
                column(width = 6, 
                       uiOutput(ns("selectSamplesUI")), 
                       uiOutput(ns("GetdFactorRef"))),
                column(width = 6, uiOutput(ns("ExpDesignTable"))),
                fluidRow(
                  
                )
              )
            #fluidRow( column(width = 6, uiOutput(ns("ValidUI"))))
            )
        ),
        fluidRow( uiOutput(outputId = ns("LoadDataUI"))),
        fluidRow( column(width = 6, actionButton(inputId = ns("loadData"),"load Data"))),
        br(),
    
        fluidRow(
          uiOutput(ns("CompletenessUI")),
          uiOutput(ns("UpsetMAE"))
        )

  )

}


LoadOmicsData <- function(input, output, session, rea.values){

    observe({

      rea.values$loadData <- FALSE
      rea.values$model    <- FALSE
      rea.values$analysis <- FALSE
      
      rea.values$contrastMat  <- NULL 
      FlomicsMultiAssay      <<- NULL
      session$userData$Design <- NULL
      
      rea.values$validate.status <- 0

    })
    #local.rea.values <- reactiveValues(Plot = FALSE)
    

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
        
        rea.values$Contrasts.Sel <- NULL 
        FlomicsMultiAssay      <<- NULL
        session$userData$Design <- NULL
        
      # rea.values$dataAnalysis   <- FALSE
      # 
      #   local.rea.values$Plot     <- FALSE
      #   rea.values$selectModel    <- FALSE
      #   rea.values$selectContrast <- FALSE
      #   rea.values$DataExplor     <- FALSE
      #   
      #   FlomicsMultiAssay <<- NULL

        # read and check design file
        ExpDesign.tbl <- read_ExpDesign(input.file = input$Experimental.Design.file$datapath)

        validate({ need(!is.null(ExpDesign.tbl), message="") })

        print("# 1- Load experimental design...")

        ####### Display tab of exp design  ########
        if(dim(ExpDesign.tbl)[2] > 5){
          ExpDesign.tbl.affich <- ExpDesign.tbl[,1:5] 
          
        }else{
          ExpDesign.tbl.affich <- ExpDesign.tbl }

        # display table
        # output$ExpDesignTable <- DT::renderDataTable({
        #   DT::datatable(ExpDesign.tbl.affich, options = list(pageLength = 5, dom = 'tp') ) })
        output$ExpDesignTable <- renderUI({
          box(width = 12, background = "light-blue",
            tags$b("The experimental design view :"),
              
            DT::renderDataTable( DT::datatable(data = ExpDesign.tbl.affich, filter = 'top', 
                                               options = list( pageLength = 5, autoWidth = TRUE, dom = 'tp' )))
                                               #caption = 'Table 1: This is a simple caption for the table.' 
            )
          })
        
        
        # output$txt <- renderText({
        #   data <- ExpDesign.tbl.affich[input$ExpDesignTable_rows_all,]
        #   print(paste0("toto",dim(data)[1])) })
        
        
        
        #### sample list  ####
        output$selectSamplesUI <- renderUI( {

          sampleList <- rownames(ExpDesign.tbl.affich)

          box(width = 12, background = "green", # Valid colors are: red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.

              # List of samples and diff modalities of factors
              tags$b("Select samples or/and modalities of factors :"),
              fluidRow(
                  column(width= ceiling(12/(length(names(ExpDesign.tbl.affich))+1)),
                    # sample list :
                    pickerInput(
                      inputId  = session$ns("selectSamples"),
                      label    = tags$span(style="color: black;","Samples"),
                      choices  = sampleList,
                      options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                      multiple = TRUE,
                      selected = sampleList)
                  ),
                  
                  # condition list
                  lapply(names(ExpDesign.tbl.affich), function(i) {
                    column(width= ceiling(12/(length(names(ExpDesign.tbl.affich))+1)),
                      # sample list :
                      pickerInput(
                        inputId  = session$ns(paste0("selectFactors.", i)),
                        label    = tags$span(style="color: black;",i),
                        choices  = levels(as.factor(ExpDesign.tbl.affich[[i]])),
                        options  = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"),
                        multiple = TRUE,
                        selected = levels(as.factor(ExpDesign.tbl.affich[[i]]))))})
              )
          )
        })

        ####### set ref and type of each factor ########
        output$GetdFactorRef <- renderUI({

          # filter samples
          ExpDesign.tbl.flt <<- ExpDesign.tbl[input$selectSamples,] 
          # filtering per conditions
          for(i in names(ExpDesign.tbl.affich)) { 
            ExpDesign.tbl.flt <<- dplyr::filter(ExpDesign.tbl.flt, get(i) %in% input[[paste0("selectFactors.", i)]])
          }
          
          box(width = 12, background = "green",

              # Construct the form to set the reference factor level
              tags$b("Select the type and the level of reference fo each design factor :"),
              fluidRow(
                lapply(names(ExpDesign.tbl.affich), function(i) {
                  column(width= round(12/length(names(ExpDesign.tbl.affich))),
                      selectInput(inputId = session$ns(paste0("dF.RefLevel.", i)), 
                                  label   = tags$span(style="color: black;",i) ,
                                  choices = unique((ExpDesign.tbl.flt[[i]]))),
                      
                      radioButtons(session$ns(paste0("dF.Type.", i)), label=NULL , choices = c("Bio","batch", "Meta"), 
                                   selected = "Bio", inline = FALSE, width = 2, choiceNames = NULL, choiceValues = NULL)
              
                      )
                  })
              )
          )
        })
        
          # output$ValidUI <- renderUI({
          #   actionButton(session$ns("ValidF"),"Valid factor set up") })

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

            uiOutput(session$ns(paste("toAddData",addDataNum + 1,sep="")))
          )
        })
    })

    # as soon as the "load" buttom has been clicked
    # => create ExpDesign object
    # => create flomicsMultiAssay object
    # => upsetR
    observeEvent(input$loadData, {

      ### load Design
      
      # reset objects and UI
      rea.values$loadData <- FALSE
      rea.values$model    <- FALSE
      rea.values$analysis <- FALSE

      rea.values$Contrasts.Sel   <- NULL 
      FlomicsMultiAssay        <<- NULL
      session$userData$Design   <- NULL
      
      # #updateTabItems(session, "tabs", selected = "importData")
      # rea.values$DataExplor     <- FALSE
      # rea.values$dataAnalysis   <- FALSE
      # 
      # local.rea.values$Plot     <- FALSE
      # rea.values$selectContrast <- FALSE

      
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
      
      # ExpDesign.tbl <- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1, sep = "\t") %>% 
      #   dplyr::mutate(dplyr::across(where(is.character), stringr::str_remove_all, pattern = fixed(" ")))
      # rownames(ExpDesign.tbl) <- stringr::str_remove_all(string = rownames(ExpDesign.tbl), pattern = " ")
      # names(ExpDesign.tbl)    <- stringr::str_remove_all(string = names(ExpDesign.tbl),    pattern = " ")
      
      ExpDesign.tbl <- ExpDesign.tbl.flt
      
      ### Get the Type and ref of the factors that the users enter in the form
      dF.Type.dFac<-vector()
      dF.List.Name<-vector()
      dF.List.ref <-vector()
      
      rea.values$validate.status <- 0
      
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
      session$userData$Design <- ExpDesign.constructor(ExpDesign = ExpDesign.tbl, refList = dF.List.ref, typeList = dF.Type.dFac)
      
      #### load omics data
      ####
      print("# 2- Load omics data...")
      inputs <- list()
      omicList <- list()
      #### list of omic data laoded from interface
      rea.values$validate.status <- 0
      dataName.vec <- c()
      for (k in 1:addDataNum){

        if(input[[paste0("omicType", k)]] != "none"){

          omicType <- input[[paste0("omicType", k)]]

          dataName <- paste0(omicType, ".", gsub("[[:space:]]", "", input[[paste0("DataName", k)]]))

          # check presence of dataname
          if(dataName == ""){
            showModal(modalDialog( title = "Error message", "Dataset names is required : dataset ", k ))
            rea.values$validate.status <- 1
          }

          # check duplicat dataset name
          if(any(duplicated(omicList[[omicType]])) == TRUE){
            showModal(modalDialog( title = "Error message", "Dataset names must be unique : dataset ", (1:addDataNum)[duplicated(dataName.vec)] ))
            rea.values$validate.status <- 1
          }

          #### read omic data file
          # => check omics data
          if(is.null(input[[paste0("data", k)]])){
            showModal(modalDialog( title = "Error message", "omics counts/abundances matrix is required : dataset ", k ))
            rea.values$validate.status <- 1
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

          validate({ need(rea.values$validate.status == 0, message="error") })
        }
      }

      FlomicsMultiAssay.try <- try_rflomics(FlomicsMultiAssay.constructor(inputs = inputs, Design=session$userData$Design, projectName=input$projectName))

      if(!is.null(FlomicsMultiAssay.try[["error"]])) {
        showModal(modalDialog( title = "Error message", as.character(FlomicsMultiAssay.try[["error"]])))
      }
      validate({ need(is.null(FlomicsMultiAssay.try[["error"]]), message="error") })


      
      FlomicsMultiAssay <<- FlomicsMultiAssay.try[["value"]]
      
      session$userData$FlomicsMultiAssay <- FlomicsMultiAssay

      # FlomicsMultiAssay[, complete.cases(FlomicsMultiAssay), ]
      # colData(FlomicsMultiAssay)[!complete.cases(FlomicsMultiAssay),]
      # intersectColumns(FlomicsMultiAssay)
     
      # => completeness
      #### check experimental design : experimental design must be a complete and balanced.
      
      print(paste0("#    => Design Completeness Check..."))

      completeCheckRes <- CheckExpDesignCompleteness(object = session$userData$Design)
      FlomicsMultiAssay@metadata[["completeCheck"]] <<- completeCheckRes
      
      if(!is.null(completeCheckRes[["error"]])){
        showModal(modalDialog(title = "Error message", completeCheckRes[["error"]]))

      }
      
      # continue only if message is true or warning
      validate({ need(is.null(completeCheckRes[["error"]]) ,message="") })
      
      # 
      rea.values$loadData <- TRUE
      rea.values$model    <- TRUE
      
    })
    
    ##########################################  
    # (3) check data
    ########################################## 

    # completeness check
    output$CompletenessUI <- renderUI({
      
      #if (local.rea.values$Plot == FALSE) return()
      if (rea.values$loadData == FALSE) return()
      
      print(paste0("#    => Completeness plot..."))
      
      box( width = 6,  status = "warning",
           
           # plot of count per condition
           renderPlot( 
             isolate({ FlomicsMultiAssay@metadata[["completeCheck"]][["plot"]] })
             ),
           hr(),
           tags$i("You **must** have a **complete design** (i.e. all possible combinations of factor's level).
                     **Balanced design** (presence of the same number of replicats for all
                     possible combinations) is not required  but advised.
                     You **must** also have at least one biological factor and 2 replicats")
      )
    })
        
    # upset of all data
    output$UpsetMAE <- renderUI({
      
     # if (local.rea.values$Plot == FALSE) return()
      if (rea.values$loadData == FALSE) return()
      if (length(unlist(FlomicsMultiAssay@metadata$omicList)) < 2) return()
      
      
      #if (length(names(FlomicsMultiAssay)) >1){
        print(paste0("#    => upset plot..."))
        
        box(width = 6, status = "warning", 
            renderPlot( isolate({ upsetSamples(FlomicsMultiAssay[unlist(FlomicsMultiAssay@metadata$omicList),]) })),
            hr(),
            tags$i("**discription**")
        )
        #}
    })
    
    return(input)
    
}


######### FUNCTION ##########

read_ExpDesign <- function(input.file){

  res <- try_rflomics(read.table(input.file, header = TRUE,row.names = 1, sep = "\t")  %>% 
         dplyr::mutate(dplyr::across(where(is.character), stringr::str_remove_all, pattern = fixed(" "))) %>% 
         dplyr::mutate(dplyr::across(.cols = where(is.character), as.factor)))
  
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

    mat <- res$value
    rownames(mat) <- stringr::str_remove_all(string = rownames(mat), pattern = " ")
    names(mat)    <- stringr::str_remove_all(string = names(mat),    pattern = " ")
    
    return(mat)
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
