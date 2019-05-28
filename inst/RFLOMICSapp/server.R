
library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

# definition of the loadData function()

  loadData <- function() {

    if(is.null(input$Count.Import.file)){
      stop("A file with count data is required as input")
    }
    else{
    count <- read.table(input$Count.Import.file$datapath,h=T,row.names = 1)
    }

    if(is.null(input$QC.Import.file)){
    QCmat <- NULL
    }
    else{
     QCmat <- read.table(input$QC.Import.file$datapath,h=T)
    }
    FE <<- FlomicsExperiment(count,QCmat)

  }

  #######
  # Part1
  #######

  # as soon as the "load" button has been clicked
  #  => the loadData function is called and the experimental design item is printed

  observeEvent(input$load, {

    loadData()

    output$ExpDesignItem <- renderMenu({
      menuItem("Experimental Design", tabName = "designExp",icon = icon('vials'))
      })

      # Construct the form to set the reference factor level
      output$GetdFactorRef <- renderUI({

                  lapply(1:FE@colDataStruc["n_dFac"], function(i) {
                    box(width=3,
                        selectInput(paste0("dF.RefLevel.dFac", i), paste0('dFac', i),
                                    choices = levels(FE@colData[,paste0('dFac', i)]),
                                    selectize=FALSE,
                                    size=5)
                    )})})

      # Construct the form to set the the type of the factor (either biological or batch)
      output$GetdFactorType <- renderUI({
                  lapply(1:FE@colDataStruc["n_dFac"], function(i) {
                    box(width=3,
                        radioButtons(paste0("dF.Type.dFac", i), label=NULL , choices = c("Bio","batch"), selected = "Bio",
                                     inline = FALSE, width = 2, choiceNames = NULL,
                                     choiceValues = NULL)
                    )})})

      # Construct the form to enter the name of the factor
      output$GetdFactorName <- renderUI({
                  lapply(1:FE@colDataStruc["n_dFac"], function(i) {
                    box(width=3,
                        textInput(paste0("dF.Name.dFac", i), label=NULL , value = paste0("dFac",i), width = NULL,
                                  placeholder = NULL)
                    )})})
      })



  # Definition of the updateDesignFactors function()
  # This function update the design ob
   updateDesignFactors <- function(){

     dF.Type.dFac<-vector()
     dF.List.Name<-vector()

     # Get the Type and the name of the factors that the users enter in the form
     for(i in 1:FE@colDataStruc["n_dFac"]){
       dF.Type.dFac[i] <- input[[paste0("dF.Type.dFac",i)]]
       dF.List.Name[i] <-input[[paste0("dF.Name.dFac",i)]]
     }

      List.Factors.new <- FE@design@List.Factors

      # Relevel the factor
       for(i in 1:length(List.Factors.new)){
         List.Factors.new[[i]] <- relevel(List.Factors.new[[i]],ref=input[[paste0("dF.RefLevel.dFac",i)]])
     }
     names(List.Factors.new) <- dF.List.Name

     FE@design@List.Factors <<- List.Factors.new
     FE@design@Factors.Type <<- dF.Type.dFac
     names(FE@colData)[1:FE@colDataStruc["n_dFac"]] <<- dF.List.Name
   }


   CheckInputFacName <- function(){
     for(i in 1:FE@colDataStruc["n_dFac"]){
       if(input[[paste0("dF.Name.dFac",i)]]==""){
         showModal(modalDialog(
           title = "Error message",
           "Empty factor are not allowed"
         ))
       }
       validate({
                 need(input[[paste0("dF.Name.dFac",i)]] != "",message="Set a name")
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

    output$ViewData <- renderMenu({

      menuItem("View Data", tabName = "Data",icon = icon('eye'), startExpanded = FALSE,
               menuSubItem("design and QC matrix",
                           tabName = "dmatrix",
                           icon = icon('table')),
               menuSubItem("count matrix",
                           tabName = "AbundanceMatrix",
                           icon = icon('table'))
      )
    })

    output$QualityCheck <- renderMenu({
      menuItem("Quality Check", tabName = "QC",icon = icon('check-square'),
               menuSubItem("Explanations",
                           tabName = "QC",
                           icon = icon('info-circle')),
               menuSubItem("Data",
                           tabName = "QCdata",
                           icon = icon('chart-bar')),
               menuSubItem("Design",
                           tabName = "QCdesign",
                           icon = icon('chart-bar'))
      )
    })

    output$colData <- renderDataTable({
      return(as.data.frame(SummarizedExperiment::colData(FE)))
    })

    output$CountMat <- renderDataTable({
      return((assay(FE)))
    })

    output$QCdesign <- renderPlot(
      mvQCdesign(FE,axis=input$nAxis)
    )

    output$QCdata <- renderPlot(
      mvQCdata(FE,axis=input$nAxis)
    )

    output$Normalize <- renderMenu({
      menuItem("Normalize",
               tabName = "Norm",badgeLabel = "new", badgeColor = "green")
    })

    # Construct the form to select the model
     output$SetModelFormula <- renderUI({
       box(width=12,selectInput( "model.formulae",
                 "Select a model formulae",
                 choices = rev(names(GetModelFormulae(names(FE@design@List.Factors)))),
                  selectize=FALSE,size=5))
   })
     output$validM <- renderUI({
     actionButton("ValidM","Valid model choice")
     })
    })


   #######
   # Part3
   #######

   # as soon as the "valid model formulae" button has been clicked
   # => The model formulae is set and the interface to select the type of contrasts appear
   # => The menu items Normalization appear


  observeEvent(input$ValidM, {
    # => Set the model formulae
    FE@design@Model.formula <- input$model.formulae

    #  => The contrasts have to be choosen
    output$SetContrasts <- renderUI({
      infoBox( "Set Contrasts", icon = icon("line-chart"),width=6)
    })
    output$Normalization<- renderMenu({
      menuItem("Normalization check", tabName = "Norm",icon = icon('check-square'))
    })
    
    #### Normalisation step
    #palette=RColorBrewer::brewer.pal(dim(FE@colData)[1], "Set1")
    
    # dans les functions
    FE <- RunNormalization(FE)

    # Boxplot check
    output$norm.boxplot <- renderPlot( boxplotQCnorm(FE) )
  
    # compute PCA
    
    ## for unnormalised data
    #pseudo <- log2(assay(FE)+1)
    #resPCA_raw <- FactoMineR::PCA(t(pseudo),ncp = 5,graph=F)
    
    ## for normalised data
    pseudo_norm <- log2(scale(assay(FE),center=FALSE,scale=FE@Normalization@Norm.factors)+1)
    FE@listPCA[["norm"]] <- FactoMineR::PCA(t(pseudo_norm),ncp = 5,graph=F)
    
    ### PCA barplot coordinates
    
    output$condColor1 <- renderUI({
      condition <- c("samples",names(FE@colData[1:FE@colDataStruc["n_dFac"]]))
      radioButtons(inputId = 'condColorSelect1', 
                   label = 'Condition :', 
                   choices = condition,
                   selected = "samples")
    })
  
    output$norm.PCAbar <- renderPlot({
      barplotPCAnorm(FE , condition=input$condColorSelect1)
    })  
 
    ### PCA point.plot coordinates
    
    # select axis to plot
    observe({
      x <- input$PC1
      # Can also set the label and select items
      choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
      updateRadioButtons(session, "PC2",
                         choices = choices[-as.numeric(x)],
                         inline = TRUE)
      })
 
    output$condColor <- renderUI({
      condition <- c("samples",names(FE@colData[1:FE@colDataStruc["n_dFac"]]))
      radioButtons(inputId = 'condColorSelect', 
                   label = 'Levels :', 
                   choices = condition,
                   selected = "samples")
      })
    
    output$norm.PCAcoord <- renderPlot({
      
      PC1.value <- as.numeric(input$PC1)
      PC2.value <- as.numeric(input$PC2)
      
      plotPCAnorm(FE, PCs=c(PC1.value, PC2.value), condition=input$condColorSelect) 
      })
    
    
    # MAplot  Normalization   
    output$NormText <- renderText("MAplot : mean expression vs log fold change")
    


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
                     nAxis = input$nAxis,
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

