
library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {


# definition of the loadData function()
# A corriger Constructor
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

    print("# 1- Load experimental design...")
    ExpDesign <<- read.table(input$Experimental.Design.file$datapath,header = TRUE,row.names = 1)

    # construct ExperimentalDesign object
    Design <<- ExperimentalDesign(ExpDesign)

  }


  loadData <- function() {

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
    FlomicsMultiAssay <<- MultiAssayExperiment(experiments = listExp,
                                               colData     = ExpDesign,
                                               sampleMap   = listToMap(listmap),
                                               metadata    = list(design = Design,
                                                                  colDataStruc = c(n_dFac = dim(ExpDesign)[2], n_qcFac = 0),
                                                                  omicList = names(listExp)))
    print("# ...FlomicsMultiAssay has been created...")
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

  #######
  # Part1
  #######

  # as soon as the "load" button has been clicked
  #  => the loadData function is called and the experimental design item is printed

  observeEvent(input$loadExpDesign, {

    loadExpDesign()
    ####### Design ########
    output$ExpDesignItem <- renderMenu({

      menuItem("Experimental Design", tabName = "designExp",icon = icon('vials'))
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
                                choices = rev(names(GetModelFormulae(Factors.Name=names(Design@List.Factors),
                                                                     Factors.Type=Design@Factors.Type))),
                                selectize=FALSE,size=5))
    })
    output$validModelFormula<- renderUI({
      actionButton("validModelFormula","Valid model choice")
    })
})

  #######
  # Part3
  #######

  # as soon as the "valid model formulae" button has been clicked
  # => The model formulae is set and the interface to select the type of contrasts appear
  # => The menu items Normalization appear


  observeEvent(input$validModelFormula, {
    # => Set the model formulae
    Design@Model.formula <<- input$model.formulae

    # => Get and Display all the contrasts
    Design@Contrasts.List <<- Contrasts.List.Bidon
    Design@Contrasts.Coeff <<- Contrasts.Coeff.bidon

    #  => The contrasts have to be choosen
    output$SetContrasts <- renderUI({
      textOutput("2 by 2 contrasts")
      box(size=3,
          checkboxGroupInput("ListOfContrasts1", "Treatment effect for each Genotype",
                             c("WT.Treated - WT.Control"="C1",
                               "dbMT.Treated - dbMT.Control"="C2",
                               "oxMT.Treated - oxMT.Control"="C3",
                               "siMT.Treated - siMT.Control"="C4")),

          checkboxGroupInput("ListOfContrasts2", "Genotype effect at each treatment",
                             c("dbMT.Control - WT.Control"="C5",
                               "siMT.Control - WT.Control"="C6",
                               "oxMT.Control - WT.Control"="C7",
                               "dbMT.Control - siMT.Control"="C8",
                               "dbMT.Control - siMT.Control"="C9",
                               "dbMT.Treated - WT.Treated"="C10",
                               "siMT.Treated - WT.Treated"="C11",
                               "oxMT.Treated - WT.Treated"="C12",
                               "dbMT.Treated - siMT.Treated"="C13",
                               "dbMT.Treated - siMT.Treated"="C14"
                             ))
      )
      #infoBox( "Set Contrasts", icon = icon("line-chart"),width=6)
    })
    output$validContrasts <- renderUI({
      actionButton("validContrasts","Valid contrast(s) choice(s)")
    })

  })


  #######
  # Part4
  #######

  # as soon as the "valid Contrasts" buttom has been clicked
  # =>
  # =>

  observeEvent(input$validContrasts, {

    Design@Contrasts.Sel <<- c(input$ListOfContrasts1,input$ListOfContrasts2)

    output$importData <- renderMenu({
      menuItem("Load Data", tabName = "importData",icon = icon('vials'))
    })

  })

    observeEvent(input$loadData, {
      ### load data
      loadData()

      ########## Item for each omics #############
      output$omics <- renderMenu({
        menu_list <- list()
        for (omic in FlomicsMultiAssay@metadata$omicList){
          menu_list <- list(
            menu_list,
            menuItem(paste0(omic, " Data Analysis"), tabName = paste0(omic,"DataAnalysis"), icon = icon('eye'), startExpanded = FALSE,
                     menuSubItem("Data Exploratory", tabName = paste0(omic,"ExploratoryQC")),
                     menuSubItem("Filter and Normalize", tabName = paste0(omic,"Normalization")),
                     menuSubItem("Differential Analysis", tabName ="DiffAnalysis")
            ))
        }
        sidebarMenu(.list = menu_list)
      })


    # library size plot
    output$LibSize <- renderPlot(
      plotLibSize(assay(FlomicsMultiAssay[["RNAseq"]])), height = 300)

    # abundance distribution
    output$CountDist <- renderPlot(
      plotDistr(assay(FlomicsMultiAssay[["RNAseq"]])), height = 300)

    FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data="RNAseq", data.norm="raw")

    ## QCdesign tabbox
    output$QCdesign <- renderPlot(
      mvQCdesign(FlomicsMultiAssay,data="RNAseq",PCA="raw", axis=5))

    ## QCdata tabbox
    output$QCdata <- renderPlot(
      mvQCdata(FlomicsMultiAssay,data="RNAseq",PCA="raw",axis=5))
   })


     # summary of loaded data

     observeEvent(input$RunFiltering, {
       ## Run Filtering
       print("# ...FilterLowAbundance...")
       FlomicsMultiAssay <<- FilterLowAbundance(FlomicsMultiAssay, data="RNAseq", input$FilterSeuil)

       output$SummaryAbundance <- renderTable({
         ## summary
         data.frame(number=c(dim(assay(FlomicsMultiAssay[["RNAseq.filtred"]])),
                             FlomicsMultiAssay@metadata$colDataStruc[1],
                             length(FlomicsMultiAssay[["RNAseq.filtred"]]@metadata$FilteredFeature)),
                    row.names=c("Features", "Samples", "Factors", "FiltredFeature"))

         #FlomicsMultiAssay[["RNAseq.filtred"]]@metadata$Normalization$coefNorm["norm.factors"]

       }, rownames=TRUE, bordered = TRUE)

     })

       observeEvent(input$RunNormalization, {

       ## Run Normalisation
       print("# ...RunNormalization...")
       FlomicsMultiAssay <<- RunNormalization(FlomicsMultiAssay, data="RNAseq.filtred", input$selectNormMethod)

       ## Boxplot check tabbox
       output$norm.boxplot <- renderPlot({
         abundanceBoxplot(FlomicsMultiAssay, dataType="RNAseq.filtred")
       })

       print("# ...runPCA...")

       ## compute PCA on filtred data

       # FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data="RNAseq.filtred", data.norm="norm")

     })

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

#     output$DiffAnalysis <- renderMenu({
#          menuItem("Differential Analysis", tabName = "DiffAnalysis", icon = icon('not-equal'), startExpanded = FALSE)
#     })
#
#     output$CoExpAnalysis <- renderMenu({
#          menuItem("Co-Expression Analysis", tabName = "CoExpAnalysis", icon = icon('project-diagram'), startExpanded = FALSE)
#     })
#

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
                                             DiffAnalysisMethod=input$AnaDiffMethod)


       output$ContrastsResults <- renderUI({

         lapply(FlomicsMultiAssay@metadata$design@Contrasts.Sel, function(i) {
           fluidRow(
            column(12,
            box(verticalLayout(renderPlot(ggplot(data=FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]])+
                                            geom_histogram(aes(x=PValue))),
                               tags$br(),
                               renderDataTable(round(FlomicsMultiAssay@ExperimentList[["RNAseq.filtred"]]@metadata[["AnaDiffDeg"]][[i]],5),
                                               options = list(rownames = TRUE,pageLength = 10))),width=12,
          # box(renderTable(FlomicsMultiAssay@metadata$design@Contrasts.List),width=10,
               solidHeader = TRUE,
               collapsible = TRUE,
               collapsed = TRUE,
               title = FlomicsMultiAssay@metadata$design@Contrasts.List[which(FlomicsMultiAssay@metadata$design@Contrasts.List$idContrast == i),]
          )
             ),
          column(1,
                 checkboxInput("checkContrasts",NULL,width=1)
          )
           )

           })})
     })

# })




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
      save(FlomicsMultiAssay,file=file.path(tempdir(), "FlomicsMultiAssay.RData"))

      # Set up parameters to pass to Rmd document
      params <- list(Count.file = input$Count.Import.file,
                     QC.file = input$QC.Import.file,
                     FEdata=file.path(tempdir(), "FlomicsMultiAssay.RData"))

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

