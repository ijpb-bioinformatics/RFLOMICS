
library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {


# definition of the loadData function()

  FlomicsSummarizedExpConstractor <- function(dataFile, qcFile){
    
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
  
  
  loadData <- function() {
    
    ### Experimental Design
    if(is.null(input$Experimental.Design.file)){
      stop("[ERROR] Experimental Design is required")
    }
    
    print("# 1- Load experimental design...")
    ExpDesign <- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1)
    
    # constract ExperimentalDesign object
    Design <- ExperimentalDesign(ExpDesign)
    
    ### list of omic data
    print("# 2- Load omic data...")
    listOmicsDataInput <- list() 
    listOmicsDataInput[["RNAseq"]]     <- list(data = input$RNAseq.Count.Import.file,      QC = input$RNAseq.QC.Import.file)
    listOmicsDataInput[["proteome"]]   <- list(data = input$prot.abundances.Import.file,   QC = input$prot.QC.Import.file)
    listOmicsDataInput[["metabolome"]] <- list(data = input$metabo.abundances.Import.file, QC = input$metabo.QC.Import.file)

    ### constract SummarisedExperiment object for each omic data
    listmap <- list() 
    listExp <- list() 
    
    for (omic in names(listOmicsDataInput)){
      
      if(! is.null(listOmicsDataInput[[omic]]$data)){
      
        print(paste0("# ...load ", omic," data..."))
        listExp[[omic]] <- FlomicsSummarizedExpConstractor(listOmicsDataInput[[omic]]$data, listOmicsDataInput[[omic]]$QC)
        listmap[[omic]] <- data.frame(primary = as.vector(listExp[[omic]]@colData$primary),
                                      colname = as.vector(listExp[[omic]]@colData$colname),
                                      stringsAsFactors = FALSE)
        }
      }
      
    # check data list
    if (is.null(listExp)){  stop("[ERROR] No data loaded !!!")  }
 
    ### constract MultiArrayExperiment object for all omic data  
    print("# ...FlomicsMultiAssay...")
    FlomicsMultiAssay <<- MultiAssayExperiment(experiments = listExp, 
                                               colData     = ExpDesign, 
                                               sampleMap   = listToMap(listmap),
                                               metadata    = list(design = Design, 
                                                                  colDataStruc = c(n_dFac = dim(ExpDesign)[2], n_qcFac = 0),
                                                                  omicList = names(listExp)))
    
  }

  #######
  # Part1
  #######

  # as soon as the "load" button has been clicked
  #  => the loadData function is called and the experimental design item is printed

  observeEvent(input$load, {

    ### load data
    loadData()
    ########## Item for each omics #############
    output$omics <- renderMenu({
      menu_list <- list()
      for (omic in FlomicsMultiAssay@metadata$omicList){
        menu_list <- list(
          menu_list,
          menuItem(paste0(omic, " Data Analysis"), tabName = paste0(omic,"DataAnalysis"), icon = icon('eye'), startExpanded = FALSE,
                menuItem("Data Exploratory", tabName = paste0(omic,"Exploratory"),startExpanded = FALSE,
                    menuSubItem("Data", tabName = paste0(omic,"ExploratoryData")),
                    menuSubItem("Bio. and tech. Variability", tabName = paste0(omic,"ExploratoryQC")))))
      }
      sidebarMenu(.list = menu_list)
    })

    ####### Design ########
    output$ExpDesignItem <- renderMenu({
    
      menuItem("Experimental Design", tabName = "designExp",icon = icon('vials'))
      })

    # Construct the form to set the reference factor level
    output$GetdFactorRef <- renderUI({

      lapply(names(FlomicsMultiAssay@colData), function(i) {
        box(width=3,
            selectInput(paste0("dF.RefLevel.", i), i,
                        choices = levels(FlomicsMultiAssay@colData[,i]),
                        selectize=FALSE,
                        size=5)
        )})})

    # Construct the form to set the the type of the factor (either biological or batch)
    output$GetdFactorType <- renderUI({
        lapply(names(FlomicsMultiAssay@colData), function(i) {
          box(width=3,
              radioButtons(paste0("dF.Type.", i), label=NULL , choices = c("Bio","batch"), selected = "Bio",
                           inline = FALSE, width = 2, choiceNames = NULL,
                           choiceValues = NULL)
          )})})

    # Construct the form to enter the name of the factor
    output$GetdFactorName <- renderUI({
                lapply(names(FlomicsMultiAssay@colData), function(i) {
                  box(width=3,
                      textInput(paste0("dF.Name.", i), label=NULL , value = i, width = NULL,
                                placeholder = NULL)
                  )})})

   })



  # Definition of the updateDesignFactors function()
  # This function update the design ob
   updateDesignFactors <- function(){

     dF.Type.dFac<-vector()
     dF.List.Name<-vector()

     # Get the Type and the name of the factors that the users enter in the form
     for(dFac in names(FlomicsMultiAssay@colData)){
       dF.Type.dFac[dFac] <- input[[paste0("dF.Type.",dFac)]]
       dF.List.Name[dFac] <- input[[paste0("dF.Name.",dFac)]]
     }

      List.Factors.new <- FlomicsMultiAssay@metadata$design@List.Factors

      # Relevel the factor
     for(dFac in names(List.Factors.new)){
         List.Factors.new[[dFac]] <- relevel(List.Factors.new[[dFac]],ref=input[[paste0("dF.RefLevel.",dFac)]])
     }
     names(List.Factors.new) <- dF.List.Name

     FlomicsMultiAssay@metadata$design@List.Factors <<- List.Factors.new
     FlomicsMultiAssay@metadata$design@Factors.Type <<- dF.Type.dFac
     names(FlomicsMultiAssay@colData)[1:FlomicsMultiAssay@metadata$colDataStruc["n_dFac"]] <<- dF.List.Name
   }


   CheckInputFacName <- function(){
     for(dFac in names(FlomicsMultiAssay@colData)){
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

   #######
   # Part2
   #######

   # as soon as the "Valid design set up" button has been clicked
   #  => The updateDesignFactors function is called
   #  => The menu items "View Data", "QualityCheck" and "Normalize" are printed
   #  => The box to choose the "model formulae" is printed


  observeEvent(input$ValidF, {

    CheckInputFacName()
    updateDesignFactors()

    # Construct the form to select the model
    output$SetModelFormula <- renderUI({
      box(width=12,selectInput( "model.formulae",
                                "Select a model formulae",
                                choices = rev(names(GetModelFormulae(Factors.Name=names(FlomicsMultiAssay@metadata$design@List.Factors),
                                                                     Factors.Type=FlomicsMultiAssay@metadata$design@Factors.Type))),
                                selectize=FALSE,size=5))
    })
    output$validM <- renderUI({
      actionButton("ValidM","Valid model choice")
    })
    

    # library size plot
    output$LibSize <- renderPlot(
      plotLibSize(assay(FlomicsMultiAssay[["RNAseq"]])), height = 300)
    
    # abundance distribution
    output$CountDist <- renderPlot(
      plotDistr(assay(FlomicsMultiAssay[["RNAseq"]])), height = 300)  
    

#   })


#   observeEvent(input$Norm, {

     
     # summary of loaded data
     output$SummaryAbundance <- renderTable({
       
       
       ## Run Filtering
       print("# ...FilterLowAbundance...")
       FlomicsMultiAssay <<- FilterLowAbundance(FlomicsMultiAssay, data="RNAseq", input$FilterSeuil)
       
       ## Run Normalisation
       print("# ...RunNormalization...")
       FlomicsMultiAssay <<- RunNormalization(FlomicsMultiAssay, data="RNAseq.filtred", input$selectNormMethod)
       
       ## summary
       data.frame(number=c(dim(assay(FlomicsMultiAssay[["RNAseq.filtred"]])), 
                           FlomicsMultiAssay@metadata$colDataStruc[1],
                           length(FlomicsMultiAssay[["RNAseq.filtred"]]@metadata$FilteredFeature)), 
                  row.names=c("Features", "Samples", "Factors", "FiltredFeature"))
       
       #FlomicsMultiAssay[["RNAseq.filtred"]]@metadata$Normalization$coefNorm["norm.factors"]
       
       }, rownames=TRUE, bordered = TRUE)

     ## Boxplot check tabbox
     output$norm.boxplot <- renderPlot({
       
       abundanceBoxplot(FlomicsMultiAssay, dataType="RNAseq.filtred")
       })

          
#     ## compute PCA on normalized data
#     pseudo_norm <- log2(scale(assay(FE),center=FALSE,scale=FE@Normalization@Norm.factors$norm.factors)+1)
#     FE@listPCA[["norm"]] <<- FactoMineR::PCA(t(pseudo_norm),ncp = 5,graph=F)
#
#     ### PCA plot
#
#     # select axis to plot
#     observe({
#       x <- input$PC1
#       # Can also set the label and select items
#       choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
#       updateRadioButtons(session, "PC2",
#                          choices = choices[-as.numeric(x)],
#                          inline  = TRUE)
#     })
#    
#     # select factors for color plot
#     output$condColor <- renderUI({
#       condition <- c("groups",names(FE@colData[1:FE@colDataStruc["n_dFac"]]))
#       radioButtons(inputId = 'condColorSelect',
#                    label = 'Levels :',
#                    choices = condition,
#                    selected = "groups")
#     })
#
#     # PCA plot
#     output$norm.PCAcoord <- renderPlot({
#
#       PC1.value <- as.numeric(input$PC1)
#       PC2.value <- as.numeric(input$PC2)
#
#       plotPCAnorm(FE, data=input$selectData2, PCs=c(PC1.value, PC2.value), condition=input$condColorSelect)
#     })
#
#     ## QCdesign tabbox
#     output$QCdesign <- renderPlot(
#       mvQCdesign(FE,data=input$selectData, axis=5))
#     
#     ## QCdata tabbox
#     output$QCdata <- renderPlot(
#       mvQCdata(FE,axis=5))
#     
#     output$DiffAnalysis <- renderMenu({
#          menuItem("Differential Analysis", tabName = "DiffAnalysis", icon = icon('not-equal'), startExpanded = FALSE)
#     })
#     
#     output$CoExpAnalysis <- renderMenu({
#          menuItem("Co-Expression Analysis", tabName = "CoExpAnalysis", icon = icon('project-diagram'), startExpanded = FALSE)
#     })
#     
 })

   #######
   # Part3
   #######

   # as soon as the "valid model formulae" button has been clicked
   # => The model formulae is set and the interface to select the type of contrasts appear
   # => The menu items Normalization appear


  observeEvent(input$ValidM, {
    # => Set the model formulae
    FlomicsMultiAssay@metadata$design@Model.formula <- input$model.formulae

    #  => The contrasts have to be choosen
    output$SetContrasts <- renderUI({
      infoBox( "Set Contrasts", icon = icon("line-chart"),width=6)
    })
    
  })


##########
# Report
##########

  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # TEST
      # save FE object in .Rdata and load it during report execution
      save(FE,file=file.path(tempdir(), "FE.RData"))

      # Set up parameters to pass to Rmd document
      params <- list(Count.file = input$Count.Import.file,
                     QC.file = input$QC.Import.file,
                     FEdata=file.path(tempdir(), "FE.RData"))

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )

})

