
library(shiny)


rm(list = ls())


shinyServer(function(input, output, session) {


    ########################################################################
    ######################### MAIN #########################################

    ##########################################
    # Part1 : Set Math model
    ##########################################
    inputExp <- callModule(ExperimentalDesign, "Exp")
   
    # display set up model Item
    observeEvent(inputExp$ValidF, {
      
      output$SetUpModelMenu <- renderMenu({
        menuSubItem("Statistical model", tabName = "SetUpModel",  selected = TRUE)
      })
      
      output$print <- renderPrint({ "toto" })
      
    })
    
    inputModel <- callModule(GLM_model, "model")
  
    observeEvent(inputModel$validContrasts, {
      
      output$importData <- renderMenu({
        menuItem("Load Data", tabName = "importData",icon = icon('download'), selected = TRUE)
      })
    })
    
    ##########################################
    # Part2 : load data
    ##########################################
    
    inputData <- callModule(LoadOmicsData, "data")

    observeEvent(inputData$loadData, {
    
          #### Item for each omics #####
          output$omics <- renderMenu({
            menu_list <- list()
            menu_list <- list(
              menu_list,
              sidebarMenu(id = "sbm",
                          
                  lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
                    
                    do.call(menuItem, c(text = paste0(omics, " Analysis"), tabName = paste0(omics, "Analysis"),
                                        
                          #lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                          
                          switch(omics ,
                                 "RNAseq"={
                                   lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                                     menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]),
                                                 tabName = paste0("RNAseqAnalysis", i), icon = icon('chart-area'), selected = FALSE)
                                   })
                                   
                                 },
                                 "proteomics"={
                                   lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                                     menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]),
                                                 tabName = paste0("ProtAnalysis", i), icon = icon('chart-area'), selected = FALSE)
                                   })
                                 },
                                 "metabolomics"={
                                   lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
                                     menuSubItem(text = paste0(FlomicsMultiAssay@metadata$omicList[[omics]][[i]]),
                                                 tabName = paste0("MetaAnalysis", i), icon = icon('chart-area'), selected = FALSE)
                                   })
                                   #menuItem(paste0(omic, " Data Exploratory"), tabName = "MetaExploratoryQC", icon = icon('chart-area'), selected = TRUE),
                                   #menuItem(paste0(omic, " Data Processing"),  tabName = "MetaProcessing",    icon = icon('chart-area'), selected = FALSE)
                                 }
                          )
                    ))
                  })
              )
            )
            sidebarMenu(.list = menu_list)
          })
    })
    


    
    
#     ##########################################
#     # Part3 : Data Exploratory
#     ##########################################
#     lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
# 
#       switch(omics ,
#              "RNAseq"={
#                  lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
#                    callModule(RNAseqDataExplorTab, paste0("RNAseq",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#                   })
#              },
#              "proteomics"={
#                  lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
#                    callModule(ProtMetaDataExplorTab, paste0("proteomics",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#                  })
#              },
#              "metabolomics"={
#                  lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
#                    callModule(ProtMetaDataExplorTab, paste0("metabolomics",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#                  })
#               }
#       )
# 
#       #lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
# 
#       #})
#     })
# 
#     ##########################################
#     # Part4 : Data processing : filtering, Normalisation...
#     ##########################################
#     lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
# 
#       #lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
#       for(i in names(FlomicsMultiAssay@metadata$omicList[[omics]])){
# 
#         callModule(RNAseqDataNormTab, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#       }
#       #})
#     })
# 
#     ##########################################
#     # Part5 : Analysi Diff
#     ##########################################
# 
#     inputDiff <- list()
#     lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
# 
#       lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
#         #callModule(DiffExpParam, i, FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#         inputDiff[[paste0(omics, i)]] <<- callModule(DiffExpAnalysis, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
# 
#       })
#     })
# 
#   ##########################################
#   # Part6 : Co-Expression Analysis
#   ##########################################
#     lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
# 
#       lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
#         observeEvent(inputDiff[[paste0(omics, i)]]$runAnaDiff, {
# 
#           callModule(CoSeqAnalysis, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#         })
#       })
#     })
#     ##########################################
#     # Part7 : Enrichment Analysis
#     ##########################################
#     lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){
# 
#       lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
# 
#         observeEvent(inputDiff[[paste0(omics, i)]]$runAnaDiff, {
# 
#           callModule(AnnotationEnrichment, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
#         })
#       })
#     })
# 
# 
# 
#   # observeEvent(input$buttonValidMerge, {
#   #
#   #   output$CoExpression <- renderMenu({
#   #     menuItem("Co-expression Analysis", tabName = "CoExpression",icon = icon('chart-area'), selected = FALSE)
#   #   })
#   #
#   #
#   #  output$Asuivre <- renderPrint({
#   #
#   #    paste0("À suivre...")
#   #  })
#   # })
# 
# # })





##########
# Report
##########

  # output$report <- downloadHandler(
  #   # For PDF output, change this to "report.pdf"
  #   filename = "report.html",
  #   content = function(file) {
  #     # Copy the report file to a temporary directory before processing it, in
  #     # case we don't have write permissions to the current working dir (which
  #     # can happen when deployed).
  # 
  #     tempReport <-  "report.Rmd" # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # 
  #     #tempReport <- file.path(tempdir(), "report.Rmd")
  #     #file.copy("report.Rmd", tempReport, overwrite = TRUE)
  # 
  #     # TEST
  #     # save FE object in .Rdata and load it during report execution
  #     save(FlomicsMultiAssay,file=file.path(tempdir(), "FlomicsMultiAssay.RData"))
  # 
  #     # Set up parameters to pass to Rmd document
  #     params <- list( FEdata = file.path(tempdir(), "FlomicsMultiAssay.RData"),
  #                     pngDir = tempdir())
  # 
  #     print(tempdir())
  #     # Knit the document, passing in the `params` list, and eval it in a
  #     # child of the global environment (this isolates the code in the document
  #     # from the code in this app).
  #     rmarkdown::render(tempReport, output_file = file,
  #                       params = params,
  #                       envir = new.env(parent = globalenv()))
  #   }
  # )

})

