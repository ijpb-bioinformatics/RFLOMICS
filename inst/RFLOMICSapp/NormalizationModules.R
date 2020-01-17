########################
# RFLOMICS MODULES
########################

# exemple
sliderTextUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  
  tagList(
    sliderInput(ns("slider"), "Slide Me", 0, 100, 1),
    textOutput(ns("number"))
  )
}

sliderText <- function(input, output, session){
  
  output$number <- renderText({ input$slider })
}



##########################################
# Part4 : Data processing : filtering, Normalisation...
##########################################

RNAseqDataNormTabUI <- function(id){
  
  #name space for id
  ns <- NS(id)
  tagList(
  
      fluidRow(
        column(4,
               fluidRow(
                 box( title = "Low Abundance Filtering",  solidHeader = TRUE, status = "warning", width = 12,
                      numericInput(inputId = ns("FilterSeuil"),
                                   label="Threshold :",
                                   value=0, 0, max=100, 1 ),
                      
                      #actionButton("RunFiltering","Run Filtering"),
                      verbatimTextOutput(ns("FilterResults"))
                 )
               ),
               fluidRow(
                 box( title = "Normalization" , solidHeader = TRUE, status = "warning", width = 12,
                      selectInput(inputId  = ns("selectNormMethod"),
                                  label    = "Method :",
                                  choices  =  list("TMM (edgeR)" = "TMM", "RLE (edgeR)" = "RLE", "upperquartile (edgeR)" = "upperquartile"),
                                  selected = "TMM"),
                      #actionButton("RunNormalization","Run Normalisation")
                 )
               ),
               actionButton(ns("NormValid"),"Validate")
        ),
        column(8,
               box(title = "Abundance distribution", solidHeader = TRUE, status = "warning", width = 14 ,  height = NULL,
                   plotOutput(ns("norm.boxplot"))
               )
        )
      ),
      fluidRow(
        
        box(title = "Principal component analysis", solidHeader = TRUE, status = "warning", width = 12 ,  height = NULL,
            
            column(width = 2,
                   fluidRow( uiOutput(ns('condColor')) ),
                   tags$br(),
                   fluidRow( uiOutput(ns('PC1axis'))),
                   fluidRow( uiOutput(ns('PC2axis'))),
                   tags$br(),
                   tags$br(),
                   fluidRow(actionButton(ns("screenshotPCA_Norm"),"Screenshot"))
            ),
            column(width = 10, plotOutput(ns("norm.PCAcoord")))
        )
      )
  )
}


RNAseqDataNormTab <- function(input, output, session, dataset){
 
  FlomicsMultiAssay.rea <<- reactive({
    #### Filter low abundance ####
    FilterSeuil <- input$FilterSeuil
    
    #### Run Normalisation ####
    print("# 8- Abundance normalization...")
    FlomicsMultiAssay <<- RunNormalization(FlomicsMultiAssay, data=paste0(dataset,".filtred"), input$selectNormMethod)
    
    #### Run PCA for filtred & normalized data ####
    FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data=paste0(dataset,".filtred"), PCA="norm")
    
    FlomicsMultiAssay
  })
  
  
  print("# 7- Low Abundance Filtering...")
  output$FilterResults <- renderPrint({
    
    FlomicsMultiAssay <<- FilterLowAbundance(FlomicsMultiAssay, data=dataset, input$FilterSeuil)
    paste0( length(FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$FilteredFeature),
            " features filtered (from ", dim(FlomicsMultiAssay[[dataset]])[1], ")")
  })
  
  ## Boxplot of distribution of normalized abundance 
  output$norm.boxplot <- renderPlot({
    FlomicsMultiAssay.rea()
    abundanceBoxplot(FlomicsMultiAssay, dataType=paste0(dataset,".filtred"), 
                     pngFile=file.path(tempdir(), paste0(dataset,"_norm.boxplot.png")))
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
    
    plotPCAnorm(FlomicsMultiAssay.rea(), data=paste0(dataset,".filtred"), PCA="norm", PCs=c(PC1.value, PC2.value), condition=input$condColorSelect, 
                file.path(tempdir(), paste0(dataset,"_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")))
  })
  
  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_Norm, {
    PC1.value <- as.numeric(input$PC1)
    PC2.value <- as.numeric(input$PC2)
    
    file.copy(file.path(tempdir(), paste0(dataset,"_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              file.path(tempdir(), paste0(dataset,"_PCAdesign_norm_PC", PC1.value,"_PC", PC2.value, "_", input$condColorSelectRaw, ".png")), 
              overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
  })
  
}



