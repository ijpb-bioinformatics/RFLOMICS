############# module 

UpdateRadioButtonsUI <- function(id){

  #name space for id
  ns <- NS(id)

  tagList(
    radioButtons(inputId  = ns("Firstaxis"),
                 label    = "Choice of PCs :",
                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                 selected = 1, inline = TRUE),

    # select PCA axis 2 for plot
    radioButtons(inputId  = ns("Secondaxis"),
                 label    = "",
                 choices  = list("PC1" = 1, "PC2" = 2, "PC3" = 3),
                 selected = 2, inline = TRUE)
  )
}


UpdateRadioButtons <- function(input, output, session){

  observeEvent(input$Firstaxis, {

    x <- input$Firstaxis
    # Can also set the label and select items
    choices=c("PC1" = 1, "PC2" = 2, "PC3" = 3)
    updateRadioButtons(session, "Secondaxis",
                       choices = choices[-as.numeric(x)],
                       inline  = TRUE)
  })

}



RadioButtonsConditionUI <- function(id){

  #name space for id
  ns <- NS(id)

  tagList(

    uiOutput(ns('condColor')),
  )
}

RadioButtonsCondition <- function(input, output, session){

  # select factors for color PCA plot
  output$condColor <- renderUI({

    condition <- c("groups",names(FlomicsMultiAssay@colData))
    radioButtons(inputId = session$ns("condColorSelect"),
                 label = 'Levels :',
                 choices = condition,
                 selected = "groups")
  })
}

###### summary of all analysed data ####
omics_data_analysis_summaryUI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    fluidRow(
      box(title = "Summary of dataset analysis", width = 12, status = "warning", solidHeader = TRUE)),
    fluidRow(
      box(title = "Dataset processing", width = 12, status = "warning", 
          column(4, plotOutput(ns("procSummary"))),
          column(8, plotOutput(ns("mofaPlot"))))),
    fluidRow(
      box(title = "Differential analysiss", width = 12, status = "warning", plotOutput(ns("DiffSummary")))),
  )
}

omics_data_analysis_summary <- function(input, output, session, rea.values){

  output$procSummary <- renderPlot({ 
    
    summaryProcess.df <- lapply(rea.values$datasetDiff, function(dataset){
      
      c(paste(dataset, "filtred", sep = "."), dim(FlomicsMultiAssay[[dataset]]) - dim(FlomicsMultiAssay[[paste(dataset, "filtred", sep = ".")]]))
    }) %>% purrr::reduce(rbind) %>% data.table()
    names(summaryProcess.df) <- c("dataset", "rm_entities", "rm_samples")
    
    ggplot(summaryProcess.df %>% melt(id="dataset"), aes(x=variable, y=dataset)) + 
      geom_tile(aes(fill=dataset),color="white", lwd =1.5, linetype=1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.ticks = element_blank(), legend.position="none") + 
      ylab("") + xlab("") + geom_text(aes(label=value)) 
    })
  
  
  output$mofaPlot <- renderPlot({ 

    MAE.filtred.Diff <- FlomicsMultiAssay[,, paste(rea.values$datasetDiff, "filtred", sep = ".")]
    MAE.MOFA <- MOFA2::create_mofa_from_MultiAssayExperiment(MAE.filtred.Diff)
    MOFA2::plot_data_overview(MAE.MOFA)
  })
  
  output$DiffSummary <- renderPlot({ 
  
    
    summaryDiff.df <- lapply( paste(rea.values$datasetDiff, "filtred", sep = "."), function(dataset){
      lapply(names(FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats), function(contrast){  
        
        c(dataset, contrast, paste0(FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats[[contrast]][["gDEup"]], "/", 
                                    FlomicsMultiAssay[[dataset]]@metadata$DiffExpAnal$stats[[contrast]][["gDEdown"]]) )
        
      }) %>% purrr::reduce(rbind)
      
    }) %>% purrr::reduce(rbind) %>% data.table()
    
    names(summaryDiff.df) <- c("dataset", "contrast", "value")
    
    
    ggplot(summaryDiff.df, aes(x=contrast, y=dataset)) + 
      geom_tile(aes(fill=dataset),color="white", lwd =1.5, linetype=1) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), axis.ticks = element_blank(), legend.position="none",
            axis.text.x = element_text(angle = 90, hjust = 1)) + 
      ylab("") + xlab("") + geom_text(aes(label=value)) + ggtitle("DE entities (up/down)")
  })
  
}


