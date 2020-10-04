
library(shiny)
# source("DataExploratoryModules.R")
# source("NormalizationModules.R")
# source("DiffExpressionModules.R")
# source("CoExpressionModules.R")

rm(list = ls())

shinyServer(function(input, output, session) {


  ############################################################
  ######################### FUNCTIONS ########################
  
  FlomicsSummarizedExpConstructor <- function(dataFile, qcFile){

    abundance <- read.table(dataFile$datapath, header = TRUE, row.names = 1)

    # RNAseq QC
    if(is.null(qcFile)){
      QCmat <- data.frame(primary = colnames(abundance),
                          colname = colnames(abundance),
                          stringsAsFactors = FALSE)
    }
    else{
      print("# ... metadata QC...")
      QCmat <- read.table(qcFile$datapath, header = TRUE)
    }
    se <- SummarizedExperiment(assays  = S4Vectors::SimpleList(abundance=as.matrix(abundance)),
                               colData = QCmat)
    #se <- SummarizedExperiment(assays  = list(counts=counts), colData = QCmat)
    return(se)
  }

  loadExpDesign <- function() {

    ### Experimental Design
    #if(is.null(input$Experimental.Design.file)){
    #  stop("[ERROR] Experimental Design is required")
    #}
    if(is.null(input$Experimental.Design.file)){
      showModal(modalDialog(
        title = "Error message",
        "Experimental Design is required"
      ))
    }
    validate({
      need(! is.null(input$Experimental.Design.file), message="Set a name")
    })
    
    
    ExpDesign <<- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1)
    
    # construct ExperimentalDesign object
    Design <<- ExperimentalDesign(ExpDesign)
    
  }

  # definition of the loadData function()
  loadData <- function() {
    
    listOmicsDataInput <<- list()
    
    #### list of omic data laoded from interface
    print("# 5- Load omic data...")
    dataName.vec <- c()
    for (k in 1:addDataNum){
      
      omicType <- input[[paste0("omicType", k)]]
      
      dataName <- paste0(omicType, ".", gsub("[[:space:]]", "", input[[paste0("DataName", k)]]))
      
      dataName.vec <- c(dataName.vec, dataName)
      
      if(omicType != "none"){
        
        # check omics data
        if(is.null(input[[paste0('data', k)]])){
          showModal(modalDialog( title = "Error message", "omics counts/abundances matrix is required : dataset ", k ))
        }
        validate({
          need(!is.null(input[[paste0('data', k)]]), message="error")
        })
        
        # check presence of dataname
        if(dataName == ""){
          showModal(modalDialog( title = "Error message", "Dataset names is required : dataset ", k ))
        }
        validate({
          need(dataName != "", message="error")
        })
        #gsub("[[:space:]]", "", x)
        
        # check duplicat dataset name
        if(any(duplicated(dataName.vec)) == TRUE){
          showModal(modalDialog( title = "Error message", "Dataset names must be unique : dataset ", (1:addDataNum)[duplicated(dataName.vec)] ))
        }
        validate({
          need(any(duplicated(dataName.vec)) == FALSE, message="error")
        })

        listOmicsDataInput[[omicType]][[dataName]] <<- list(data = input[[paste0("data", k)]], 
                                                            QC   = input[[paste0("metadataQC", k)]],
                                                            dataName = dataName,
                                                            omicType = omicType,
                                                            order    = k)
      }
    } 
    
  }
  
  
  FlomicsMultiAssayExperimentConstructor <- function(){
    
    listmap <- list()
    listExp <- list()
    omicList<- list()
    
    for (omicType in names(listOmicsDataInput)){
    
      for (dataName in names(listOmicsDataInput[[omicType]])){
      
        ### build SummarisedExperiment object for each omic data
        if(! is.null(listOmicsDataInput[[omicType]][[dataName]]$data)){
          
          colnames <- c(names(omicList[[omicType]]), listOmicsDataInput[[omicType]][[dataName]]$order)
          omicList[[omicType]] <- c(omicList[[omicType]], dataName)
          names(omicList[[omicType]]) <- colnames

          print(paste0("# ...load ", dataName," data..."))
          listExp[[dataName]] <- FlomicsSummarizedExpConstructor(listOmicsDataInput[[omicType]][[dataName]]$data, 
                                                                 listOmicsDataInput[[omicType]][[dataName]]$QC)
          listmap[[dataName]] <- data.frame(primary = as.vector(listExp[[dataName]]@colData$primary),
                                            colname = as.vector(listExp[[dataName]]@colData$colname),
                                            stringsAsFactors = FALSE)
        }
      }
    }
    
    # check data list
    if (is.null(listExp)){  stop("[ERROR] No data loaded !!!")  }

    ### constract MultiArrayExperiment object for all omic data
    print("# 6- Build FlomicsMultiAssay object...")
    FlomicsMultiAssay <<- MultiAssayExperiment(experiments = listExp,
                                               colData     = ExpDesign,
                                               sampleMap   = listToMap(listmap),
                                               metadata    = list(design = Design,
                                                                  colDataStruc = c(n_dFac = dim(ExpDesign)[2], n_qcFac = 0),
                                                                  omicList = omicList))
  }

  # Definition of the updateDesignFactors function()
  # This function update the design ob
  updateDesignFactors <- function(){

    dF.Type.dFac<-vector()
    dF.List.Name<-vector()
    
    old.names <- names(Design@ExpDesign)
    current.names <- names(Design@List.Factors)

    # Get the Type and the name of the factors that the users enter in the form
    for(dFac in old.names){
      
      dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
      dF.List.Name[dFac] <- input[[paste0("dF.Name.",dFac)]]
    }

    List.Factors.new <- Design@List.Factors

    # Relevel the factor
    for(dFac.n in 1:length(old.names)){
      
      tmp <- List.Factors.new[[current.names[dFac.n]]] %>% as.factor()
      List.Factors.new[[current.names[dFac.n]]] <- relevel(tmp, ref=input[[paste0("dF.RefLevel.",old.names[dFac.n])]])
    }
    names(List.Factors.new) <- dF.List.Name
    names(dF.Type.dFac)     <- dF.List.Name

    Design@List.Factors <<- List.Factors.new
    Design@Factors.Type <<- dF.Type.dFac
    #names(Design@List.Factors)[1:length(Design@List.Factors)] <<- dF.List.Name
  }


  CheckInputFacName <- function(){
    # check if empty
    for(dFac in names(Design@ExpDesign)){
      if(input[[paste0("dF.Name.",dFac)]]==""){
        showModal(modalDialog(
          title = "Error message",
          "Empty factor are not allowed"
        ))
      }
      validate({
        need(input[[paste0("dF.Name.",dFac)]] != "",message="Set a name")
      })
    }
  }

  
  ########################################################################
  ######################### MAIN #########################################
  
  ##########################################
  # Part2 : Define the Experimental design:
  #         -> load experimental plan
  #         -> the level of ref for each factor
  #         -> the formulae
  #         -> the model
  #         -> the contrasts
  ##########################################

  # as soon as the "load" button has been clicked
  #  => the loadExpDesign function is called and the experimental design item is printed
  observeEvent(input$loadExpDesign, {
    
    print("# 1- Load experimental design...")
    
    #shinyjs::reset("SetUpModel")
    #resetForm(SetUpModel)
    loadExpDesign()

    # display desgin table
    output$ExpDesignTable <- renderUI({

        box(width = 8, status = "warning",
            DT::renderDataTable( DT::datatable(ExpDesign) )
          )
      })
  
  ####### Set up Design model ########
  output$GetdFactorRef <- renderUI({

    box(status = "warning", width = 12, height = NULL,
        
       # Construct the form to set the reference factor level
       h4("Select the level of reference fo each design factor"),
       fluidRow(
         lapply(names(Design@List.Factors), function(i) {
           box(width=3,
               selectInput(paste0("dF.RefLevel.", i), i,
                           choices = levels(Design@List.Factors[[i]]),
                           selectize=FALSE,
                           size=5))})
       ),
       
       # Construct the form to set the type of the factor (either biological or batch)
       h4("Select the type of the design factor"),
       fluidRow( 
         lapply(names(Design@List.Factors), function(i) {
           box(width=3,
               radioButtons(paste0("dF.Type.", i), label=NULL , choices = c("Bio","batch"), selected = "Bio",
                            inline = FALSE, width = 2, choiceNames = NULL,
                            choiceValues = NULL))})
         
       ),
       
        # Construct the form to enter the name of the factor
        h4("Enter a name for each design factor"),
        fluidRow( 
          lapply(names(Design@ExpDesign), function(i) {
            box(width=3,
                textInput(paste0("dF.Name.", i), label=NULL , value = i, width = NULL,
                          placeholder = NULL))})
          
        ),

       actionButton("ValidF","Valid factor set up")
     )
    
    })
  })
  
  # as soon as the "Valid factor set up" button has been clicked
  #  => The upvdateDesignFactors function is called
  #  => The interface to select the model formulae appear
  observeEvent(input$ValidF, {
    print("# 2- Set design model...")
    
    # 
    CheckInputFacName()
    updateDesignFactors()

    #### check experimental design : experimental design must be a complete and balanced.
    completeCheckRes <- CheckExpDesignCompleteness(Design)
    
    output$Completeness <- renderUI({
      
      column(width= 12, 
              box( status = "warning", width = 12, 
      
                  # print message
                  renderText( completeCheckRes[["message"]][2] ),
                  
                  # plot of count per condition
                  renderPlot( plotExperimentalDesign(completeCheckRes[["count"]] ))
                  
              )
        )
    })
    
    
    ## error/warning message
    if(completeCheckRes[["message"]][1] == "false"){
      showModal(modalDialog(
        title = "Error message",
        completeCheckRes[["message"]][2]
      ))
    }
    
    if(completeCheckRes[["message"]][1] == "warning"){
      showModal(modalDialog(
        title = "Warning message",
        completeCheckRes[["message"]][2]
      ))
    }
    
  
    ## continue only if message is true or warning
    validate({
      need(completeCheckRes[["message"]][1] != "false" ,message="ok")
    })
    
    
    # display set up model Item
    output$SetUpModel <- renderMenu({
      menuSubItem("Design matrix", tabName = "SetUpModel",  selected = TRUE)
      
    })
    
    # Construct the form to select the model
    output$SetModelFormula <- renderUI({
      box(status = "warning", width = 12,
          h4("Select a model formulae"),
          selectInput( "model.formulae", "",
                        choices = rev(names(GetModelFormulae(Factors.Name=names(Design@List.Factors),
                                                             Factors.Type=Design@Factors.Type))),
                        selectize=FALSE,size=5),
          actionButton("validModelFormula","Valid model choice")
          )
      })
    })

  # as soon as the "valid model formulae" button has been clicked
  # => The model formulae is set and the interface to select the contrasts appear
  observeEvent(input$validModelFormula, {
    
    print("# 3- Choice of statistical model...")

    # => Set the model formulae
    Design@Model.formula <<- input$model.formulae
    print(paste0("#    model :", Design@Model.formula))
    
    # => get list of expression contrast (hypothesis)
    Design <<- getExpressionContrast(Design)
    
    
    # => Set Model design matrix
    # => Get and Display all the contrasts
    
    
    #  => The contrasts have to be choosen
    output$SetContrasts <- renderUI({
      #textOutput("2 by 2 contrasts")
      box(width=12, status = "warning", size=3,
          
      lapply(names(Design@Contrasts.List), function(contrastType) {
        
            vect        <- as.vector(Design@Contrasts.List[[contrastType]]$contrast)
            names(vect) <- as.vector(Design@Contrasts.List[[contrastType]]$contrastName)

            #checkboxGroupInput("ListOfContrasts1", paste0(i, " effect"), vect)
            box(
              checkboxGroupInput(paste0("ContrastType",contrastType), paste0("Contrast type : ", contrastType), vect)
            )
        }),

        column(width=4, actionButton("validContrasts","Valid contrast(s) choice(s)")))
      })
    })
  
  
  # as soon as the "valid Contrasts" buttom has been clicked
  # => The selected contrasts are saved
  # => The load data item appear
  observeEvent(input$validContrasts, {

    # get list of selected contrast data frames with expression, name and type 
    contrastList <- list()
    contrastList <- lapply(names(Design@Contrasts.List), function(contrastType) {
        
        #contrastList<-input[[paste0("ContrastType",contrastType)]]
        tmp <- Design@Contrasts.List[[contrastType]] %>% 
          dplyr::filter(contrast %in% input[[paste0("ContrastType",contrastType)]]) %>% 
          dplyr::select(contrast, contrastName, type, groupComparison)
        return(tmp)
    })
    Design@Contrasts.Sel <<- contrastList %>% purrr::reduce(rbind)
    
    # check if user has selected the contrasts to test
    if(dim(Design@Contrasts.Sel)[1] == 0){
      
      showModal(modalDialog(
        title = "Error message",
        "Please select the hypotheses to test."
      ))
    }
    
    ## continue only if message is true
    validate({
      need(dim(Design@Contrasts.Sel)[1] != 0, message="ok")
    })
    
    
    # define all the coefficients of selected contrasts and return a contrast matrix with contrast sample name and associated coefficients
    Design <<- getContrastMatrix(Design)
    
    output$contrastVec <- renderPrint({
      
      Design@Contrasts.Coeff
    })
   
    
    
    # 
    output$importData <- renderMenu({
      menuItem("Load Data", tabName = "importData",icon = icon('download'), selected = TRUE)
      })
    })
 
  
  
  ##########################################
  # Part2 : load data
  ##########################################
  
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
                   selectInput(inputId=paste0('omicType', addDataNum), label='Omics', 
                               choices = c("None"="none", "RNAseq"="RNAseq", 
                                           "Proteomics"="proteomics", "Metabolomics"="metabolomics"),
                               selected = "none")
                   ),
            column(4,
                   # matrix count/abundance input
                   fileInput(inputId=paste0("data", addDataNum), "omic count/abundance (Ex.)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                   ),
            column(4,
                   # metadata/QC bioinfo
                   fileInput(inputId=paste0("metadataQC", addDataNum), "QC or metadata (Ex.)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                   ),
            column(2,
                   # dataset Name
                   textInput(inputId=paste0("DataName", addDataNum), label="Dataset name", value=paste0("set", as.character(addDataNum)))
                   )
            ),
        uiOutput(paste("toAddData",addDataNum + 1,sep=""))
        )
      })
    })
  
  

  # as soon as the "load data" buttom has been clicked
  # => selected input(s) are loaded with loadData() function
  # => create tmp dir for png report
  observeEvent(input$loadData, {
    
    ### load data
    loadData()
    
    ### load data
    FlomicsMultiAssayExperimentConstructor()
    
    
    ##### data list
    #choiceList <- FlomicsMultiAssay@metadata$omicList %>% purrr::reduce(c)
    #output$dataList <- renderUI({
    #  selectInput(inputId='datalist', label='Dataset list :', 
    #              choices = choiceList,
    #              selected = as.character(choiceList[1]))
    #  })
    
    
    
    #### Item for each omics ####
    #observeEvent(input$datalist, {
    
    output$omics <- renderMenu({
      menu_list <- list()
      menu_list <- list(
        menu_list,
        sidebarMenu(id = "sbm",
                    
            lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
              
              do.call(menuItem, c(text = paste0(omics, " Analysis"), tabName = paste0(omics, "Analysis"), 
                                  
                  #lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                    
                    switch(omics ,
                           "RNAseq"={
                             lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                               menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]), 
                                        tabName = paste0("RNAseqAnalysis", i), icon = icon('chart-area'), selected = FALSE)
                             })
                            
                           },
                           "proteomics"={
                             lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                               menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]), 
                                         tabName = paste0("ProtAnalysis", i), icon = icon('chart-area'), selected = FALSE)
                             })
                           },
                           "metabolomics"={
                             lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                               menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]), 
                                         tabName = paste0("MetaAnalysis", i), icon = icon('chart-area'), selected = FALSE)
                             })
                               #menuItem(paste0(omic, " Data Exploratory"), tabName = "MetaExploratoryQC", icon = icon('chart-area'), selected = TRUE),
                               #menuItem(paste0(omic, " Data Processing"),  tabName = "MetaProcessing",    icon = icon('chart-area'), selected = FALSE)
                           }
                    )
                  #})
              ))
            })
        )
      )
      sidebarMenu(.list = menu_list)
    })
    
    
    ##########################################
    # Part3 : Data Exploratory
    ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
      
      switch(omics ,
             "RNAseq"={
                 lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                    
                   callModule(RNAseqDataExplorTab, paste0("RNAseq",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
                  })
             },
             "proteomics"={
                 lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                   
                   callModule(ProtMetaDataExplorTab, paste0("proteomics",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
                 })             
             },
             "metabolomics"={
                 lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                   
                   callModule(ProtMetaDataExplorTab, paste0("metabolomics",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
                 })
              }
      )
      
      #lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
        
        
      #})
    })
    
    ##########################################
    # Part4 : Data processing : filtering, Normalisation...
    ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
      
      #lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
      for(i in names(FlomicsMultiAssay@metadata$omicList[[omics]])){
        
        callModule(RNAseqDataNormTab, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
      }
      #})
    })

    ##########################################
    # Part5 : Analysi Diff
    ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
      
      lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
        
        #callModule(DiffExpParam, i, FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
        
        callModule(DiffExpAnalysis, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
        
      })
    })
  
    })

  
  ##########################################
  # Part6 : Co-Expression Analysis
  ##########################################
  
  observeEvent(input$buttonValidMerge, {
  
    output$CoExpression <- renderMenu({
      menuItem("Co-expression Analysis", tabName = "CoExpression",icon = icon('chart-area'), selected = FALSE)
    })
    
    
   output$Asuivre <- renderPrint({
  
     paste0("À suivre...")
   })
  })

# })




##########
# Report
##########

  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      
      tempReport <-  "report.Rmd" # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      #tempReport <- file.path(tempdir(), "report.Rmd")
      #file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # TEST
      # save FE object in .Rdata and load it during report execution
      save(FlomicsMultiAssay,file=file.path(tempdir(), "FlomicsMultiAssay.RData"))

      # Set up parameters to pass to Rmd document
      params <- list( FEdata = file.path(tempdir(), "FlomicsMultiAssay.RData"),
                      pngDir = tempdir())

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )

})

