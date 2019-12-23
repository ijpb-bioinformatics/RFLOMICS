
library(shiny)


shinyServer(function(input, output, session) {




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
    if(is.null(input$Experimental.Design.file)){
      stop("[ERROR] Experimental Design is required")
    }


    ExpDesign <<- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1)

    # construct ExperimentalDesign object
    Design <<- ExperimentalDesign(ExpDesign)

  }

  # definition of the loadData function()

  loadData <- function() {

    ### list of omic data
    print("# 5- Load omic data...")
    listOmicsDataInput <- list()
    listOmicsDataInput[["RNAseq"]]       <- list(data = input$RNAseq.Count.Import.file,      QC = input$RNAseq.QC.Import.file)
    listOmicsDataInput[["proteomics"]]   <- list(data = input$prot.abundances.Import.file,   QC = input$prot.QC.Import.file)
    listOmicsDataInput[["metabolomics"]] <- list(data = input$metabo.abundances.Import.file, QC = input$metabo.QC.Import.file)

    ### constract SummarisedExperiment object for each omic data
    listmap <- list()
    listExp <- list()

    for (omic in names(listOmicsDataInput)){

      if(! is.null(listOmicsDataInput[[omic]]$data)){

        print(paste0("# ...load ", omic," data..."))
        listExp[[omic]] <- FlomicsSummarizedExpConstructor(listOmicsDataInput[[omic]]$data, listOmicsDataInput[[omic]]$QC)
        listmap[[omic]] <- data.frame(primary = as.vector(listExp[[omic]]@colData$primary),
                                      colname = as.vector(listExp[[omic]]@colData$colname),
                                      stringsAsFactors = FALSE)
        }
      }

    # check data list
    if (is.null(listExp)){  stop("[ERROR] No data loaded !!!")  }

    ### constract MultiArrayExperiment object for all omic data
    print("# ...FlomicsMultiAssay...")
    omicList <- as.list(names(listExp))
    names(omicList) <- names(listExp)
    FlomicsMultiAssay <<- MultiAssayExperiment(experiments = listExp,
                                               colData     = ExpDesign,
                                               sampleMap   = listToMap(listmap),
                                               metadata    = list(design = Design,
                                                                  colDataStruc = c(n_dFac = dim(ExpDesign)[2], n_qcFac = 0),
                                                                  omicList = omicList))
  }

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

  ##########################################
  # Part1 : load Experimental design
  #########################################

  # as soon as the "load" button has been clicked
  #  => the loadExpDesign function is called and the experimental design item is printed

  observeEvent(input$loadExpDesign, {

    print("# 1- Load experimental design...")
    loadExpDesign()


    output$ExpDesignTable <- renderUI({

        box(width = 8, status = "warning",
            DT::renderDataTable( DT::datatable(ExpDesign) )
        )
    })

    ####### Set up Design model ########

    output$SetUpModel <- renderMenu({

      menuSubItem("Design matrix", tabName = "SetUpModel",  selected = TRUE)
    })

    # Construct the form to set the reference factor level
    output$GetdFactorRef <- renderUI({

      lapply(names(Design@List.Factors), function(i) {
        box(width=3,
            selectInput(paste0("dF.RefLevel.", i), i,
                        choices = levels(Design@List.Factors[[i]]),
                        selectize=FALSE,
                        size=5)
        )})})

    # Construct the form to set the type of the factor (either biological or batch)
    output$GetdFactorType <- renderUI({
      lapply(names(Design@List.Factors), function(i) {
        box(width=3,
            radioButtons(paste0("dF.Type.", i), label=NULL , choices = c("Bio","batch"), selected = "Bio",
                         inline = FALSE, width = 2, choiceNames = NULL,
                         choiceValues = NULL)
        )})})

    # Construct the form to enter the name of the factor
    output$GetdFactorName <- renderUI({
      lapply(names(Design@List.Factors), function(i) {
        box(width=3,
            textInput(paste0("dF.Name.", i), label=NULL , value = i, width = NULL,
                      placeholder = NULL)
        )})})

  })



  # Definition of the updateDesignFactors function()
  # This function update the design ob


  ##########################################
  # Part2 : Define the Experimental design:
  #         -> the level of ref for each factor
  #         -> the formulae
  #         -> the model
  #         -> the contrasts
  #########################################


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
    print(paste0("model :", Design@Model.formula))
    Design <<- SetModelMatrix(Design)


    # => Get and Display all the contrasts
    # =>
    #tmp <- apply (combn(colnames(A),2), 2, function(x) {
    #
    #  paste(x, collapse=" - ")
    #})




    #Design@Contrasts.List <<- Contrasts.List.Bidon
    #Design@Contrasts.Coeff <<- Contrasts.Coeff.bidon

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


  #######
  # Part3
  #
  # Load data
  #
  #######

  # as soon as the "valid Contrasts" buttom has been clicked
  # => The selected contrasts are saved
  # => The load data item appear

  observeEvent(input$validContrasts, {

    Design@Contrasts.Sel <<- c(input$ListOfContrasts1)
    
    #step_tag <<- step_tag+1

    output$importData <- renderMenu({
      menuItem("Load Data", tabName = "importData",icon = icon('download'), selected = TRUE)
    })

    ## activate next item
    #newtab <- switch(input$StateSave,
    #                 "SetUpModel" = "importData",
    #                 "importData" = "SetUpModel"
    #)
    #updateTabItems(session, "StateSave", newtab)

  })

    observeEvent(input$loadData, {
      ### load data
      loadData()

      ########## Item for each omics #############
      output$omics <- renderMenu({
        menu_list <- list()
        for (omic in FlomicsMultiAssay@metadata$omicList){
          dir.create(file.path(tmpDir, omic))
          dir.create(file.path(tmpDir, omic, "images"))
          dir.create(file.path(tmpDir, omic, "tables"))
          dir.create(file.path(tmpDir, omic, "RData"))
          dir.create(file.path(tmpDir, omic, "tmp"))
          
          menu_list <- list(
            menu_list,
            menuItem(paste0(omic, " Data Analysis"),   tabName = paste0(omic,"DataAnalysis"), icon = icon('chart-area'), startExpanded = FALSE,
                     menuSubItem("Data Exploratory",    tabName = paste0(omic,"ExploratoryQC")),
                     menuSubItem("Filter and Normalize", tabName = paste0(omic,"Normalization")),
                     menuSubItem("Differential Analysis", tabName = "DiffAnalysis"),
                     menuSubItem("Co-expression Analysis", tabName = "CoExpression")
            ))
        }
        sidebarMenu(.list = menu_list)
      })
      

    # library size plot
    output$LibSize <- renderPlot(height = 300, {
      
      plotLibSize(abundances=assay(FlomicsMultiAssay[["RNAseq"]]), pngDir=file.path(tmpDir, "RNAseq/images"))
      })
      
    # abundance distribution
    output$CountDist <- renderPlot(height = 300, {
      
      plotDistr(assay(FlomicsMultiAssay[["RNAseq"]]), pngDir=file.path(tmpDir, "RNAseq/images"))
      })

    # run PCA
    FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data="RNAseq", PCA="raw")

    ### PCA plot
    output$PCA1axisRaw <- renderUI({
      radioButtons(inputId  = "PC1raw",
                   label    = "Choice of PCs :",
                   choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                   selected = 1, inline = TRUE)
    })

    output$PCA2axisRaw <- renderUI({
      radioButtons(inputId  = "PC2raw",
                   label    = "",
                   choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                   selected = 2, inline = TRUE)
    })

    # select axis to plot
    observeEvent(input$PC1raw, {
      x <- input$PC1raw
      # Can also set the label and select items
      choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
      updateRadioButtons(session, "PC2raw",
                         choices = choices[-as.numeric(x)],
                         inline  = TRUE)
    })


    # select axis to plot
    #observe({
    #  x <- input$PC1raw
    #  # Can also set the label and select items
    #  choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
    #  updateRadioButtons(session, "PC2raw",
    #                     choices = choices[-as.numeric(x)],
    #                     inline  = TRUE)
    #})

    # select factors for color plot
    output$condColorRaw <- renderUI({
      condition <- c("groups",names(FlomicsMultiAssay@colData))
      radioButtons(inputId = 'condColorSelectRaw',
                   label = 'Levels :',
                   choices = condition,
                   selected = "groups")
    })

    # PCA plot
    output$QCdesignPCARaw <- renderPlot({

      PC1.value <- as.numeric(input$PC1raw)
      PC2.value <- as.numeric(input$PC2raw)   
      plotPCAnorm(FlomicsMultiAssay, data="RNAseq", PCA="raw", PCs=c(PC1.value, PC2.value), condition=input$condColorSelectRaw, pngDir=file.path(tmpDir, "RNAseq/tmp"))
      
      })
    
    # save current PCA plot with fixed PC
    observeEvent(input$screenshotPCA_QC, {
    
      PC1.value <- as.numeric(input$PC1raw)
      PC2.value <- as.numeric(input$PC2raw) 

      filename = paste0("PCAdesign_" , "RNAseq" , "_PC", PC1.value, "-PC", PC2.value, "_", input$condColorSelectRaw, ".png")
      
      file.copy(file.path(tmpDir, "RNAseq","tmp", filename), file.path(tmpDir, "RNAseq","images", filename), 
                overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
      
      })

    ## QCdesign distance plot
    output$QCdesignPCA <- renderPlot({
      
      mvQCdesign(FlomicsMultiAssay,data="RNAseq",PCA="raw", axis=5, pngDir=file.path(tmpDir, "RNAseq/images")) })

    ## QCdata tabbox
    output$QCdata <- renderPlot({
      mvQCdata(FlomicsMultiAssay,data="RNAseq",PCA="raw",axis=5, pngDir=file.path(tmpDir, "RNAseq/images")) })

  })

     # summary of loaded data

     observeEvent(input$RunFiltering, {
       ## Run Filtering
       print("# ...FilterLowAbundance...")
       FlomicsMultiAssay <<- FilterLowAbundance(FlomicsMultiAssay, data="RNAseq", input$FilterSeuil)

       output$FilterResults <- renderPrint({

         paste0( length(FlomicsMultiAssay[["RNAseq.filtred"]]@metadata$FilteredFeature),
                 " features filtered (from ", dim(FlomicsMultiAssay[["RNAseq"]])[1], ")")
       })

     })

       observeEvent(input$RunNormalization, {

       ## Run Normalisation
       print("# 6- RNAseq normalization...")
       FlomicsMultiAssay <<- RunNormalization(FlomicsMultiAssay, data="RNAseq.filtred", input$selectNormMethod)

       ## Boxplot check tabbox
       output$norm.boxplot <- renderPlot({
         abundanceBoxplot(FlomicsMultiAssay, dataType="RNAseq.filtred", pngDir=file.path(tmpDir, "RNAseq/images"))
       })

       print("# ...runPCA...")

       ## compute PCA on filtred data
       FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data="RNAseq.filtred", PCA="norm")

       ### PCA plot
       output$PC1axis <- renderUI({
           radioButtons(inputId  = "PC1",
                        label    = "Choice of PCs :",
                        choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                        selected = 1, inline = TRUE)
        })

        output$PC2axis <- renderUI({
           radioButtons(inputId  = "PC2",
                        label    = "",
                        choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                        selected = 2, inline = TRUE)
        })

       # select axis to plot
        observeEvent(input$PC1, {
         x <- input$PC1
         # Can also set the label and select items
         choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
         updateRadioButtons(session, "PC2",
                            choices = choices[-as.numeric(x)],
                            inline  = TRUE)
       })

       # select factors for color plot
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

         plotPCAnorm(FlomicsMultiAssay, data="RNAseq.filtred", PCA="norm", PCs=c(PC1.value, PC2.value), condition=input$condColorSelect, pngDir=file.path(tmpDir, "RNAseq/tmp"))
         
         })
       
       observeEvent(input$screenshotPCA_Norm, {
         PC1.value <- as.numeric(input$PC1)
         PC2.value <- as.numeric(input$PC2)
         
         filename = paste0("PCAdesign_" , "RNAseq.filtred" , "_PC", PC1.value, "-PC", PC2.value, "_", input$condColorSelect, ".png")
         
         file.copy(file.path(tmpDir, "RNAseq","tmp", filename), file.path(tmpDir, "RNAseq","images", filename), 
                   overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
       })
       

  })
  
     observeEvent(input$runAnaDiff, {

      # Run the differential analysis for each contrast set
      #   -> return a dynamic user interface with a collapsible box for each contrast
      #         - Pvalue graph
      #         - MAplot
      #         - Table of the DE genes

       print("# ...RunDiffAnalysis...") 

       FlomicsMultiAssay <<- RunDiffAnalysis(FlomicsMultiAssay,
                                             data="RNAseq.filtred",
                                             FDR=input$FDRSeuil ,
                                             DiffAnalysisMethod=input$AnaDiffMethod,
                                             clustermq=input$clustermq)


       output$ContrastsResults <- renderUI({

         vect <- as.vector(FlomicsMultiAssay@metadata$design@Contrasts.List$hypoth)
         names(vect) <- as.vector(FlomicsMultiAssay@metadata$design@Contrasts.List$idContrast)
         
         lapply(FlomicsMultiAssay@metadata$design@Contrasts.Sel, function(i) {
           
           resTable <- FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]]
           
           fluidRow(
            column(10,
              box(width=12, solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "warning", 
                  title = paste0(i, " : ", vect[i]),
                  
                  verticalLayout(
                    
                    ### pvalue plot ###
                    renderPlot({
                      
                      pvalue.plot(data=FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]], 
                                  tag=gsub(" ", "", vect[i]), pngDir=file.path(tmpDir, "RNAseq", "images"))
                      }),
                    #renderPlot({
                    #
                    #  data <- FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]]
                    #  # significant is coded as TRUE, not sig as FALSE
                    #  data$sig <- as.factor(abs(data$logFC) > input$logFCSeuil & data$FDR < input$FDRSeuil)
                    #  #convert FDR to -log(FDR)
                    #  data$negLogFDR <- -log10(data$FDR)
                    #
                    #  ggplot(data,aes(x=logCPM, y=logFC, color=sig)) +
                    #   geom_point() +
                    #   coord_cartesian() +
                    #   ylab("log2 FC") +
                    #   xlab("log2 CPM")
                    #
                    #}),
                    tags$br(),
                    DT::renderDataTable({
                      
                      DT::datatable(round(resTable[resTable$FDR <= input$FDRSeuil,],5), options = list(rownames = FALSE, pageLength = 10))
                      })
                    
                    # Button
                    #downloadButton("downloadData", "Download")
                    
                    # Downloadable csv of selected dataset ----
                    #output$downloadData <- downloadHandler(
                    #  filename = function() {
                    #    paste(input$dataset, ".csv", sep = "")
                    #  },
                    #  content = function(file) {
                    #    write.csv(datasetInput(), file, row.names = FALSE)
                    #  }
                    #)
                    
                    #write.csv(x = round(resTable[resTable$FDR <= input$FDRSeuil,],5), 
                    #          file = paste0(tmpDir, "/tables/DEG_",i, "_",  gsub(" ", "", vect[i]), "_" ,input$FDRSeuil, "_", input$AnaDiffMethod, ".png" ))
                    
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

         mat2venn[[i]][["features"]] <-  row.names(FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]])
         mat2venn[[i]][[i]] <- rep(1, dim(FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]])[1])
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
                              inline = FALSE, width = 2, choiceNames = NULL, choiceValues = "NULL"),
                 actionButton("buttonValidMerge","Valid")
              )
             )
           )
         )
       })
     })

     ########################################
     ##
     ########################################

     observeEvent(input$buttonValidMerge, {

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
      #save(FlomicsMultiAssay,file=file.path(tmpDir, "FlomicsMultiAssay.RData"))
      
      # Set up parameters to pass to Rmd document
      params <- list(   tag = step_tag,
                     pngDir = tmpDir,
                     FEdata = file.path(tempdir(), "FlomicsMultiAssay.RData"))

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )

})

