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
               actionButton(ns("normUpdate"),"Update")
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
                   fluidRow(
                     RadioButtonsConditionUI(ns("normData"))),
                   tags$br(),
                   UpdateRadioButtonsUI(ns("normData")),
                   tags$br(),
                   tags$br(),
                   fluidRow(actionButton(ns("screenshotPCA_Norm"),"Screenshot"))
            ),
            column(width = 10, plotOutput(ns("norm.PCAcoord")))
        )
      )
  )
}

RunFilterNormPCAfunction <- function(dataset, FilterSeuil, NormMethod){
  #### Filter low abundance ####
  print("# 7- Low Abundance Filtering...")
  FlomicsMultiAssay <- FilterLowAbundance(FlomicsMultiAssay, data=dataset, FilterSeuil)

  #### Run Normalisation ####
  print("# 8- Abundance normalization...")
  FlomicsMultiAssay <- RunNormalization(FlomicsMultiAssay, data=paste0(dataset,".filtred"), NormMethod)

  #### Run PCA for filtred & normalized data ####
  FlomicsMultiAssay <- RunPCA(FlomicsMultiAssay, data=paste0(dataset,".filtred"), PCA="norm")

  return(FlomicsMultiAssay)
}

RNAseqDataNormTab <- function(input, output, session, dataset){
  
  # FlomicsMultiAssay.rea <<- reactive({
  #   #### Filter low abundance ####
  #   print("# 7- Low Abundance Filtering...")
  #   FlomicsMultiAssay <<- FilterLowAbundance(FlomicsMultiAssay, data=dataset, FilterSeuil)
  # 
  #   #### Run Normalisation ####
  #   print("# 8- Abundance normalization...")
  #   FlomicsMultiAssay <<- RunNormalization(FlomicsMultiAssay, data=paste0(dataset,".filtred"), NormMethod)
  # 
  #   #### Run PCA for filtred & normalized data ####
  #   FlomicsMultiAssay <<- RunPCA(FlomicsMultiAssay, data=paste0(dataset,".filtred"), PCA="norm")
  # 
  #   FlomicsMultiAssay
  # })
  
  FlomicsMultiAssay <<- RunFilterNormPCAfunction (dataset, FilterSeuil=0 , NormMethod="TMM")
  
  ## Boxplot of distribution of normalized abundance 
  output$norm.boxplot <- renderPlot({
    abundanceBoxplot(FlomicsMultiAssay, dataType=paste0(dataset,".filtred"), 
                     pngFile=file.path(tempdir(), paste0(dataset,"_norm.boxplot.png")))
  })  
  

  #### PCA analysis ####
  # select PCA axis for plot
  # update/adapt PCA axis
  callModule(UpdateRadioButtons, "normData")
  
  # select factors for color PCA plot
  callModule(RadioButtonsCondition, "normData")
  
  
  # PCA plot
  output$norm.PCAcoord <- renderPlot({
    
    PC1.value <- as.numeric(input$`normData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`normData-Secondaxis`[1])   
    condGroup <- input$`normData-condColorSelect`[1]
    
    plotPCAnorm(FlomicsMultiAssay, data=paste0(dataset,".filtred"), PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup, 
                pngFile=file.path(tempdir(), paste0(dataset,".filtred","_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")))
  })
    
  observeEvent(input$normUpdate, {
    
    FilterSeuil <- input$FilterSeuil
    NormMethod  <- input$selectNormMethod
    
    
    FlomicsMultiAssay <<- RunFilterNormPCAfunction (dataset, FilterSeuil , NormMethod)
    
    output$FilterResults <- renderPrint({
  
      paste0( length(FlomicsMultiAssay[[paste0(dataset,".filtred")]]@metadata$FilteredFeature),
              " features filtered (from ", dim(FlomicsMultiAssay[[dataset]])[1], ")")
    })
    
    ## Boxplot of distribution of normalized abundance 
    output$norm.boxplot <- renderPlot({
      abundanceBoxplot(FlomicsMultiAssay, dataType=paste0(dataset,".filtred"), 
                       pngFile=file.path(tempdir(), paste0(dataset,"_norm.boxplot.png")))
    })
    
    
    # PCA plot
    output$norm.PCAcoord <- renderPlot({
      
      PC1.value <- as.numeric(input$`normData-Firstaxis`[1])
      PC2.value <- as.numeric(input$`normData-Secondaxis`[1])   
      condGroup <- input$`normData-condColorSelect`[1]
      
      plotPCAnorm(FlomicsMultiAssay, data=paste0(dataset,".filtred"), PCA="norm", PCs=c(PC1.value, PC2.value), condition=condGroup, 
                  pngFile=file.path(tempdir(), paste0(dataset,".filtred","_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")))
    })
  })
  
  
  # save current PCA plot with fixed axix & color
  ## screenShot
  observeEvent(input$screenshotPCA_Norm, {
    PC1.value <- as.numeric(input$`normData-Firstaxis`[1])
    PC2.value <- as.numeric(input$`normData-Secondaxis`[1])
    condGroup <- input$`normData-condColorSelect`[1]
    
    file.copy(file.path(tempdir(), paste0(dataset,".filtred","_PCAdesign_norm_tmp_PC", PC1.value,"_PC", PC2.value, "_", condGroup, ".png")), 
              file.path(tempdir(), paste0(dataset,".filtred","_PCAdesign_norm_PC",     PC1.value,"_PC", PC2.value, "_", condGroup, ".png")), 
              overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
    })
  
}



