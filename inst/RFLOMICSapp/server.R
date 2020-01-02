
library(shiny)
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
                                                            omicType = omicType)
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
          
          omicList[[omicType]] <- c(omicList[[omicType]], dataName)

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

    # Get the Type and the name of the factors that the users enter in the form
    for(dFac in names(Design@List.Factors)){
      dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
      dF.List.Name[dFac] <- input[[paste0("dF.Name.",dFac)]]

    }

    List.Factors.new <- Design@List.Factors

    # Relevel the factor
    for(dFac in names(List.Factors.new)){
      List.Factors.new[[dFac]] <- relevel(List.Factors.new[[dFac]],ref=input[[paste0("dF.RefLevel.",dFac)]])
    }
    names(List.Factors.new) <- dF.List.Name

    Design@List.Factors <<- List.Factors.new
    Design@Factors.Type <<- dF.Type.dFac
    names(Design@List.Factors)[1:length(Design@List.Factors)] <<- dF.List.Name
  }


  CheckInputFacName <- function(){
    for(dFac in names(Design@List.Factors)){
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
    
    loadExpDesign()

    # display desgin table
    output$ExpDesignTable <- renderUI({

        box(width = 8, status = "warning",
            DT::renderDataTable( DT::datatable(ExpDesign) )
            )
      })
    
    # display set up model Item
    output$SetUpModel <- renderMenu({
      menuSubItem("Design matrix", tabName = "SetUpModel",  selected = TRUE)
      
      })
    
    })
  
  
  ####### Set up Design model ########
  
  # Construct the form to set the reference factor level
  output$GetdFactorRef <- renderUI({

    lapply(names(Design@List.Factors), function(i) {
      box(width=3,
          selectInput(paste0("dF.RefLevel.", i), i,
                      choices = levels(Design@List.Factors[[i]]),
                      selectize=FALSE,
                      size=5)
          )
      })
    })

  # Construct the form to set the type of the factor (either biological or batch)
  output$GetdFactorType <- renderUI({
    lapply(names(Design@List.Factors), function(i) {
      box(width=3,
          radioButtons(paste0("dF.Type.", i), label=NULL , choices = c("Bio","batch"), selected = "Bio",
                       inline = FALSE, width = 2, choiceNames = NULL,
                       choiceValues = NULL)
          )
      })
    })

  # Construct the form to enter the name of the factor
  output$GetdFactorName <- renderUI({
    lapply(names(Design@List.Factors), function(i) {
      box(width=3,
          textInput(paste0("dF.Name.", i), label=NULL , value = i, width = NULL,
                    placeholder = NULL)
          )
      })
    })

  # as soon as the "Valid factor set up" button has been clicked
  #  => The upvdateDesignFactors function is called
  #  => The interface to select the model formulae appear
  observeEvent(input$ValidF, {
    print("# 2- Set design model...")
    CheckInputFacName()
    updateDesignFactors()

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

    # => Set Model design matrix
    # => Get and Display all the contrasts
    print(paste0("#    model :", Design@Model.formula))
    Design <<- SetModelMatrix(Design)


    #  => The contrasts have to be choosen
    output$SetContrasts <- renderUI({
      #textOutput("2 by 2 contrasts")
      box(width=12, status = "warning", size=3,
      lapply(unique(Design@Contrasts.List$factors), function(i) {
            vect <- as.vector(filter(Design@Contrasts.List, factors==i)[["idContrast"]])
            names(vect) <- as.vector(filter(Design@Contrasts.List, factors==i)[["hypoth"]])

            checkboxGroupInput("ListOfContrasts1", paste0(i, " effect"), vect)
        }),

        column(width=4, actionButton("validContrasts","Valid contrast(s) choice(s)")))
      })
    })
  
  
  # as soon as the "valid Contrasts" buttom has been clicked
  # => The selected contrasts are saved
  # => The load data item appear
  observeEvent(input$validContrasts, {

    Design@Contrasts.Sel <<- c(input$ListOfContrasts1)
    
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
    choiceList <- FlomicsMultiAssay@metadata$omicList %>% purrr::reduce(c)
    output$dataList <- renderUI({
      selectInput(inputId='datalist', label='Dataset list :', 
                  choices = choiceList,
                  selected = choiceList[1])
      })
    
    })

  #### selected data ####
  datasetInput <- reactive({
    
    input$datalist
    })
    
  #### Item for each omics ####
  observeEvent(input$datalist, {
    
    omic <- na.exclude(str_extract(datasetInput(), SupportedOmics))[[1]]
    
    output$omics <- renderMenu({
      menu_list <- list()

      switch(omic ,
             "RNAseq"={
                 menu_list <- list(
                   menu_list,
                   menuItem(paste0(omic, " Data Exploratory"),     tabName = "RNAseqExploratoryQC", icon = icon('chart-area'), selected = TRUE),
                   menuItem(paste0(omic, " Filter and Normalize"), tabName = "RNAseqNormalization", icon = icon('chart-area'), selected = FALSE)
                   )
               },
             "proteomics"={
                 menu_list <- list(
                   menu_list,
                   menuItem(paste0(omic, " Data Exploratory"), tabName = "ProtExploratoryQC", icon = icon('chart-area'), selected = TRUE),
                   menuItem(paste0(omic, " Data Processing"),  tabName = "ProtProcessing", icon = icon('chart-area'), selected = FALSE)
                   )
               },
             "metabolomics"={
                 menu_list <- list(
                   menu_list,
                   menuItem(paste0(omic, " Data Exploratory"), tabName = "MetaExploratoryQC", icon = icon('chart-area'), selected = TRUE),
                   menuItem(paste0(omic, " Data Processing"),  tabName = "MetaProcessing",    icon = icon('chart-area'), selected = FALSE)
                   )
               }
             )
      sidebarMenu(.list = menu_list)

    })
    
    # run PCA with raw data
    FlomicsMultiAssay <<-  RunPCA(FlomicsMultiAssay, data=datasetInput(), PCA="raw")
    
    # initialize ui contrastResults for new dataset
    output$ContrastsResults <- renderUI({})
    output$ResultsMerge     <- renderUI({})
    #removeUI(selector = paste0("#ContrastsResults"))
  })
    
  ##########################################
  # Part3 : Data Exploratory
  ##########################################
  
  #### library size plot #### 
  output$LibSize <- renderPlot(height = 300, {
    
    plotLibSize(abundances=assay(FlomicsMultiAssay[[datasetInput()]]), dataName=datasetInput(), pngFile=file.path(tempdir(), paste0(datasetInput(),"_LibSize.png")))
    })
  
  #### abundance distribution #### 
  output$CountDist <- renderPlot(height = 300, {
    
    plotDistr(abundances=assay(FlomicsMultiAssay[[datasetInput()]]), dataName=datasetInput(), pngFile=file.path(tempdir(), paste0(datasetInput(),"_CountDist.png")))
    })
  
  #### PCA analysis ####
  # select PCA axis 1 for plot
  output$PCA1axisRaw <- renderUI({
    
    radioButtons(inputId  = "PC1raw",
                 label    = "Choice of PCs :",
                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                 selected = 1, inline = TRUE)
    })
  
  # select PCA axis 2 for plot
  output$PCA2axisRaw <- renderUI({
    
    radioButtons(inputId  = "PC2raw",
                 label    = "",
                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                 selected = 2, inline = TRUE)
    })
  
  # update/adapt PCA axis
  observeEvent(input$PC1raw, {
    
    x <- input$PC1raw
    # Can also set the label and select items
    choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
    updateRadioButtons(session, "PC2raw",
                       choices = choices[-as.numeric(x)],
                       inline  = TRUE)
    })

  
  # select factors for color PCA plot
  output$condColorRaw <- renderUI({
    
    condition <- c("groups",names(FlomicsMultiAssay@colData))
    radioButtons(inputId = 'condColorSelectRaw',
                 label = 'Levels :',
                 choices = condition,
                 selected = "groups")
    })
  
  # run PCA plot
  output$QCdesignPCARaw <- renderPlot({
    
    PC1.value <- as.numeric(input$PC1raw)
    PC2.value <- as.numeric(input$PC2raw)   
    plotPCAnorm(FlomicsMultiAssay, data=datasetInput(), PCA="raw", PCs=c(PC1.value, PC2.value), 
                condition=input$condColorSelectRaw, 
                pngFile=file.path(tempdir(), paste0(datasetInput(),"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")))
    })
  

  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_QC, {
    
    PC1.value <- as.numeric(input$PC1raw)
    PC2.value <- as.numeric(input$PC2raw) 
    
    file.copy(file.path(tempdir(), paste0(datasetInput(),"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              file.path(tempdir(), paste0(datasetInput(),"_PCAdesign_raw_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
    })
  
  #### PCA analysis QCdesign ####
  output$QCdesignPCA <- renderPlot({
    
    mvQCdesign(FlomicsMultiAssay,data=datasetInput(),PCA="raw", axis=5, 
               pngFile=file.path(tempdir(), paste0(datasetInput(),"_PCAdesignCoordRaw.png"))) 
    })

  #### PCA analysis QCdata ####
  output$QCdata <- renderPlot({
    
    mvQCdata(FlomicsMultiAssay,data=datasetInput(),PCA="raw",axis=5, 
             pngFile=file.path(tempdir(), paste0(datasetInput(),"_PCAmetaCorrRaw.png"))) 
    })
  
 
  ##########################################
  # Part4 : Data processing : filtering, Normalisation...
  ##########################################
  
  FlomicsMultiAssay.rea <<- reactive({
    #### Filter low abundance ####
    FilterSeuil <- input$FilterSeuil
    
    #### Run Normalisation ####
    print("# 8- Abundance normalization...")
    FlomicsMultiAssay <<- RunNormalization(FlomicsMultiAssay, data=paste0(datasetInput(),".filtred"), input$selectNormMethod)
    
    #### Run PCA for filtred & normalized data ####
    FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data=paste0(datasetInput(),".filtred"), PCA="norm")
    
    FlomicsMultiAssay
    })
  
  
  print("# 7- Low Abundance Filtering...")
  output$FilterResults <- renderPrint({
    
    FlomicsMultiAssay <<- FilterLowAbundance(FlomicsMultiAssay, data=datasetInput(), input$FilterSeuil)
    paste0( length(FlomicsMultiAssay[[paste0(datasetInput(),".filtred")]]@metadata$FilteredFeature),
            " features filtered (from ", dim(FlomicsMultiAssay[[datasetInput()]])[1], ")")
    })

  ## Boxplot of distribution of normalized abundance 
  output$norm.boxplot <- renderPlot({
    FlomicsMultiAssay.rea()
    abundanceBoxplot(FlomicsMultiAssay, dataType=paste0(datasetInput(),".filtred"), 
                     pngFile=file.path(tempdir(), paste0(datasetInput(),"_norm.boxplot.png")))
  })
  
    
  #### PCA analysis ####
  # select PCA axis 1 for plot
  output$PC1axis <- renderUI({
     radioButtons(inputId  = "PC1",
                  label    = "Choice of PCs :",
                  choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                  selected = 1, inline = TRUE)
   })
  
  # select PCA axis 2 for plot
  output$PC2axis <- renderUI({
     radioButtons(inputId  = "PC2",
                  label    = "",
                  choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                  selected = 2, inline = TRUE)
   })

  # update/adapt PCA axis
  observeEvent(input$PC1, {
   x <- input$PC1
   # Can also set the label and select items
   choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
   updateRadioButtons(session, "PC2",
                      choices = choices[-as.numeric(x)],
                      inline  = TRUE)
   })

  # select factors for color PCA plot
  output$condColor <- renderUI({
   condition <- c("groups",names(FlomicsMultiAssay@colData))
   radioButtons(inputId = 'condColorSelect',
                label = 'Levels :',
                choices = condition,
                selected = "groups")
   })
  
  # PCA plot
  output$norm.PCAcoord <- renderPlot({
  
    PC1.value <- as.numeric(input$PC1)
    PC2.value <- as.numeric(input$PC2)
    
    plotPCAnorm(FlomicsMultiAssay.rea(), data=paste0(datasetInput(),".filtred"), PCA="norm", PCs=c(PC1.value, PC2.value), condition=input$condColorSelect, 
               file.path(tempdir(), paste0(datasetInput(),"_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")))
    })
  
  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_Norm, {
    PC1.value <- as.numeric(input$PC1)
    PC2.value <- as.numeric(input$PC2)
    
    file.copy(file.path(tempdir(), paste0(datasetInput(),"_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              file.path(tempdir(), paste0(datasetInput(),"_PCAdesign_norm_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
    })
    
  
  ##########################################
  # Part5 : Analysi Diff
  ##########################################
  
  observeEvent(input$NormValid, {
    
    output$DiffAnalysis <- renderMenu({
      menuItem("Differential Analysis", tabName = "DiffAnalysis",icon = icon('chart-area'), selected = FALSE)
    })
    output$ContrastsResults <- renderUI({})
    
  })
  
  output$DiffParam <- renderUI({
    
        box(title = datasetInput(), width = 12, status = "warning",
            column(5,
                   selectInput("AnaDiffMethod", label = "Method :",
                               choices = list("glmfit (edgeR)"="edgeRglmfit",
                                              "limma" = "limma"),
                               selected = "edgeRglmfit")
            ),
            column(3,
                   numericInput(inputId = "FDRSeuil",
                                label="FDR :",
                                value=0.05, 0, max=1, 0.01)
            ),
            column(3,
                   selectInput(inputId = "clustermq",
                               label="send job to cluster",
                               choices = list("no"=FALSE,"genotoul"=TRUE))
            ),
            column(5,
                   actionButton("runAnaDiff","Run the differential analysis")
            )
        )
      
    })
  
  
  # Run the differential analysis for each contrast set
  #   -> return a dynamic user interface with a collapsible box for each contrast
  #         - Pvalue graph
  #         - MAplot
  #         - Table of the DE genes 
  observeEvent(input$runAnaDiff, {
    
    print("# 9- Diff Analysis...") 
  
    # run diff analysis with select method
    FlomicsMultiAssay <<- RunDiffAnalysis(FlomicsMultiAssay,
                                          data=paste0(datasetInput(),".filtred"),
                                          FDR =input$FDRSeuil ,
                                          DiffAnalysisMethod=input$AnaDiffMethod,
                                          clustermq=input$clustermq)
 
    output$ContrastsResults <- renderUI({
    
        vect <- as.vector(FlomicsMultiAssay@metadata$design@Contrasts.List$hypoth)
        names(vect) <- as.vector(FlomicsMultiAssay@metadata$design@Contrasts.List$idContrast)
       
        lapply(FlomicsMultiAssay@metadata$design@Contrasts.Sel, function(i) {
         
          resTable <- FlomicsMultiAssay@ExperimentList[[paste0(datasetInput(),".filtred")]]@metadata[["AnaDiffDeg"]][[i]]
          
          fluidRow(
            column(10,
              box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", title = paste0(i, " : ", vect[i]),
                  
                  verticalLayout(
                    
                    ### pvalue plot ###
                    renderPlot({
                      
                      pvalue.plot(data    =FlomicsMultiAssay@ExperimentList[[paste0(datasetInput(),".filtred")]]@metadata[["AnaDiffDeg"]][[i]],
                                  contrast=as.vector(filter(FlomicsMultiAssay@metadata$design@Contrasts.List , idContrast == i)$hypoth),
                                  pngFile =file.path(tempdir(), paste0(datasetInput(), "_PvalueDistribution_", gsub(" ", "", vect[i]), ".png")))
                      }),
                      tags$br(),
                      DT::renderDataTable({
                        
                        DT::datatable(round(resTable[resTable$FDR <= input$FDRSeuil,],5), options = list(rownames = FALSE, pageLength = 10))
                        })
                  
                      # Button
                      #downloadButton("downloadData", "Download")
                    )
                  )
              ),
            column(2, checkboxInput(inputId = "checkContrasts", label = "validate" ,value = TRUE , width=1))
            )
          })
      })


    # merge diff results
    mat2venn <- list()
    for(i in FlomicsMultiAssay@metadata$design@Contrasts.Sel) {
    
      mat2venn[[i]][["features"]] <-  row.names(FlomicsMultiAssay@ExperimentList[[paste0(datasetInput(),".filtred")]]@metadata[["AnaDiffDeg"]][[i]])
      mat2venn[[i]][[i]] <- rep(1, dim(FlomicsMultiAssay@ExperimentList[[paste0(datasetInput(),".filtred")]]@metadata[["AnaDiffDeg"]][[i]])[1])
      mat2venn[[i]] <- tbl_df(mat2venn[[i]])
      }
      
    mat2venn.df <- mat2venn %>% purrr::reduce(dplyr::full_join, by="features")
    
    mat2venn.df[is.na(mat2venn.df)] <- 0
    
    output$ResultsMerge <- renderUI({
      fluidRow(
       column(10,
         box(width=12, solidHeader = TRUE, title = "Combine results", collapsible = TRUE, collapsed = TRUE, status = "warning",
          column(width = 7,
             renderPlot({
               title(main = "snp")
               ven::venn(mat2venn.df[,-1] , ilab=TRUE, zcolor = "style")
               })
             ),
          column(width = 5,
             radioButtons("choise", label="" , choices = c("union","intersect"), selected = "union",
                          inline = FALSE, width = 2, choiceNames = NULL, choiceValues = NULL),
             actionButton("buttonValidMerge","Valid")
             )
          )
         )
       )
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

