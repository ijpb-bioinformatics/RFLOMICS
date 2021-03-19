
##########################################
# Part3 : Data Exploratory
##########################################


##########
# RNAseq
##########

RNAseqDataExplorTabUI <- function(id){

  #name space for id
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "Raw Data Summary", solidHeader = TRUE, status = "warning", width = 12, height = NULL,
          # library size plot
          column(6, plotOutput(ns("LibSize"),   height = "400%")),
          # count distribution plot
          column(6, plotOutput(ns("CountDist"), height = "400%"))
      )
    ),
    fluidRow(
      box(title = "Exploratory of Biological and Technical variability", solidHeader = TRUE, width = 12, status = "warning",
          tabBox( id = "ExplorAnalysisQC", width = 12,

                  tabPanel("Principal component analysis (1/2)",
                           tags$br(),
                           tags$br(),
                           column(width = 2,

                                  fluidRow(
                                    RadioButtonsConditionUI(ns("rawData"))),
                                  tags$br(),
                                  UpdateRadioButtonsUI(ns("rawData")),
                                  tags$br(),
                                  tags$br(),
                                  fluidRow(actionButton(ns("screenshotPCA_QC"),"Screenshot"))
                           ),
                           column(width = 10,  plotOutput(ns("QCdesignPCARaw")))

                  ),
                  tabPanel("Principal component analysis (2/2)", plotOutput(ns("QCdesignPCA")))#,
                  #tabPanel("Quality check for technical issues", plotOutput(ns("QCdata")))
          )
      )
    )
  )
}



RNAseqDataExplorTab <- function(input, output, session, dataset){


  #### library size plot ####
  output$LibSize <- renderPlot(height = 300, {

    plotLibSize(abundances=assay(FlomicsMultiAssay@ExperimentList[[dataset]]))
  })

  #### abundance distribution ####
  output$CountDist <- renderPlot(height = 300, {

    plotDistr(abundances=assay(FlomicsMultiAssay@ExperimentList[[dataset]]))
  })

  #### PCA analysis ####
  # select PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "rawData")

  # select factors for color PCA plot
  callModule(RadioButtonsCondition, "rawData")

  # run PCA plot
  output$QCdesignPCARaw <- renderPlot({
    FlomicsMultiAssay@ExperimentList[[dataset]] <<-  RunPCA(FlomicsMultiAssay@ExperimentList[[dataset]])
    PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    condGroup <- input$`rawData-condColorSelect`[1]

    plotPCA(FlomicsMultiAssay@ExperimentList[[dataset]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
  })


  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_QC, {

    PC1.value <- as.numeric(input$`rawData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawData-Secondaxis`[1])
    condGroup <- input$`rawData-condColorSelect`[1]

    # file.copy(file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
    #           file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
    #           overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
  })

  #### PCA analysis QCdesign ####
  output$QCdesignPCA <- renderPlot({

    mvQCdesign(FlomicsMultiAssay,data=dataset,PCA="raw", axis=5,
               pngFile=file.path(tempdir(), paste0(dataset,"_PCAdesignCoordRaw.png")))
  })

  #### PCA analysis QCdata ####
  output$QCdata <- renderPlot({

    mvQCdata(FlomicsMultiAssay,data=dataset,PCA="raw",axis=5,
             pngFile=file.path(tempdir(), paste0(dataset,"_PCAmetaCorrRaw.png")))
  })

}



##########
# Proteomic
##########


ProtMetaDataExplorTabUI <- function(id){

  #name space for id
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "Raw Data Summary", solidHeader = TRUE, status = "warning", width = 12, height = NULL,
          column(6,
                 # Nombre de prot
                 strong("Number of Proteins:"),
                 verbatimTextOutput(ns("NbProt")),
                 # Nombre de NA
                 br(),
                 strong("Number of Proteins with at least 1 NA:"),
                 verbatimTextOutput(ns("NbNA")),
                 br(),
                 # Data transformation ?
                 radioButtons(
                   inputId  ="dataTransform",
                   "Which transformation did you apply to the data ?",
                   c("none" = "none",
                     "log2" = "log2",
                     "log10" = "log10")
                 )
                 ),
          # count distribution plot
          column(6, plotOutput(ns("CountDistbis"), height = "400%"))
      )
    )
    ,
    fluidRow(
      box(title = "Exploratory of Biological and Technical variability", solidHeader = TRUE, width = 12, status = "warning",
          tabBox( id = "ExplorAnalysisQC", width = 12,

                  tabPanel("Principal component analysis (1/2)",
                           tags$br(),
                           tags$br(),
                           column(width = 2,

                                  fluidRow(
                                    RadioButtonsConditionUI(ns("rawDatabis"))),
                                  tags$br(),
                                  UpdateRadioButtonsUI(ns("rawDatabis")),
                                  tags$br(),
                                  tags$br(),
                                  fluidRow(actionButton(ns("screenshotPCA_QCbis"),"Screenshot"))
                           ),
                           column(width = 10,  plotOutput(ns("QCdesignPCARawbis")))

                  ),
                  tabPanel("Principal component analysis (2/2)", plotOutput(ns("QCdesignPCAbis"))),
                  tabPanel("Quality check for technical issues", plotOutput(ns("QCdatabis")))
          )
      )
    )
  )
}



ProtMetaDataExplorTab <- function(input, output, session, dataset){

  #### abundance distribution ####
  output$CountDistbis <- renderPlot(height = 300, {

    #plotDistr(abundances=assay(FlomicsMultiAssay[[dataset]]), dataName=dataset, pngFile=file.path(tempdir(), paste0(dataset,"_CountDist.png")))
    plotDistr(abundances=assay(FlomicsMultiAssay@ExperimentList[[dataset]]))
  })


  #### Nombre de Prot, Nombre de NA
  NbProt <- dim(assays(FlomicsMultiAssay@ExperimentList[[dataset]])$abundance)[1]
  output$NbProt <- renderPrint({ NbProt })

  #### Nombre de proteine avec au moins une valeur manquante
  ## Nb NA
  NbProtWoutNA <- dim(na.omit(assays(FlomicsMultiAssay@ExperimentList[[dataset]])$abundance))[1]
  output$NbNA <- renderPrint({ NbProt - NbProtWoutNA})

  #### Data Transformation:
  observeEvent(input$dataTransform,{


  })

  #### PCA analysis ####
  # select PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "rawDatabis")

  # select factors for color PCA plot
  callModule(RadioButtonsCondition, "rawDatabis")

  # run PCA plot
  output$QCdesignPCARawbis <- renderPlot({
    FlomicsMultiAssay@ExperimentList[[dataset]] <<-  RunPCA(FlomicsMultiAssay@ExperimentList[[dataset]])
    PC1.value <- as.numeric(input$`rawDatabis-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawDatabis-Secondaxis`[1])
    condGroup <- input$`rawDatabis-condColorSelect`[1]

    #plotPCA(FlomicsMultiAssay@ExperimentList[[dataset]], data=dataset, PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup,
    #            pngFile=file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")))
    plotPCA(FlomicsMultiAssay@ExperimentList[[dataset]], PCA="raw", PCs=c(PC1.value, PC2.value), condition=condGroup)
  })


  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_QCbis, {

    PC1.value <- as.numeric(input$`rawDatabis-Firstaxis`[1])
    PC2.value <- as.numeric(input$`rawDatabis-Secondaxis`[1])
    condGroup <- input$`rawDatabis-condColorSelect`[1]

    # file.copy(file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
    #           file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")),
    #           overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
  })

  #### PCA analysis QCdesign ####
  output$QCdesignPCAbis <- renderPlot({

    mvQCdesign(FlomicsMultiAssay,data=dataset,PCA="raw", axis=5,
               pngFile=file.path(tempdir(), paste0(dataset,"_PCAdesignCoordRaw.png")))
  })

  #### PCA analysis QCdata ####
  output$QCdatabis <- renderPlot({

    mvQCdata(FlomicsMultiAssay,data=dataset,PCA="raw",axis=5,
             pngFile=file.path(tempdir(), paste0(dataset,"_PCAmetaCorrRaw.png")))
  })

}



##########
# Metabolomic
##########


