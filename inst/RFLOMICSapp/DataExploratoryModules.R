
##########################################
# Part3 : Data Exploratory
##########################################

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
                                  
                                  fluidRow(uiOutput(ns('condColorRaw'))),
                                  tags$br(),
                                  fluidRow(uiOutput(ns('PCA1axisRaw'))),
                                  fluidRow(uiOutput(ns('PCA2axisRaw'))),
                                  tags$br(),
                                  tags$br(),
                                  fluidRow(actionButton(ns("screenshotPCA_QC"),"Screenshot"))
                           ),
                           column(width = 10,  plotOutput(ns("QCdesignPCARaw")))
                           
                  ),
                  tabPanel("Principal component analysis (2/2)", plotOutput(ns("QCdesignPCA"))),
                  tabPanel("Quality check for technical issues", plotOutput(ns("QCdata")))
          )
      )
    )
  )
}



RNAseqDataExplorTab <- function(input, output, session, dataset){
  
  
  #### library size plot #### 
  output$LibSize <- renderPlot(height = 300, {
    
    plotLibSize(abundances=assay(FlomicsMultiAssay[[dataset]]), dataName=dataset, pngFile=file.path(tempdir(), paste0(dataset,"_LibSize.png")))
  })
  
  #### abundance distribution #### 
  output$CountDist <- renderPlot(height = 300, {
    
    plotDistr(abundances=assay(FlomicsMultiAssay[[dataset]]), dataName=dataset, pngFile=file.path(tempdir(), paste0(dataset,"_CountDist.png")))
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
    FlomicsMultiAssay <<-  RunPCA(FlomicsMultiAssay, data=dataset, PCA="raw")
    PC1.value <- as.numeric(input$PC1raw)
    PC2.value <- as.numeric(input$PC2raw)   
    plotPCAnorm(FlomicsMultiAssay, data=dataset, PCA="raw", PCs=c(PC1.value, PC2.value), 
                condition=input$condColorSelectRaw, 
                pngFile=file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")))
  })
  
  
  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_QC, {
    
    PC1.value <- as.numeric(input$PC1raw)
    PC2.value <- as.numeric(input$PC2raw) 
    
    file.copy(file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              file.path(tempdir(), paste0(dataset,"_PCAdesign_raw_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
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