
library(shiny)
library(shinydashboard)
library(shinyFiles)

rm(list = ls())


shinyServer(function(input, output, session) {


    ########################################################################
    ######################### MAIN #########################################

    ##########################################
    # Part0 : presentation page
    ##########################################
    # à faire
    output$print <- renderPrint({ "HELLO" })

    # # dir
    # shinyDirChoose(input = input, id = 'dir0', roots = c(home = '~'))
    # output$filepaths <- renderPrint({parseDirPath(roots = c(home = '~'), selection = input$dir0)})



    ##########################################
    # Part1 : Set GLM model
    ##########################################
    # load design
    # set reference
    # set type of factor (bio/batch)
    # check design (complete and balanced)
    inputExp <- callModule(ExperimentalDesign, "Exp")

    # display set up model Item
    # if no error message
    observeEvent(inputExp$ValidF, {

      #continue only if message is true or warning
      validate({
        need(validate.status == 0 ,message="set design step failed")
      })

      output$SetUpModelMenu <- renderMenu({
        menuSubItem("Statistical model", tabName = "SetUpModel",  selected = TRUE)
      })
    })

    # set GLM model
    # and select list of contrast to test
    inputModel <- callModule(GLM_model, "model")

    # display load data Item
    # if no error message
    observeEvent(inputModel$validContrasts, {

      #continue only if message is true
      validate({
        need(validate.status == 0 ,message="select model and contrast step failed")
      })

      output$importData <- renderMenu({
        menuItem("Load Data", tabName = "importData",icon = icon('download'), selected = TRUE)
      })
    })


    ##########################################
    # Part2 : load data
    ##########################################

    # load omics data
    inputData <- callModule(LoadOmicsData, "data")

    # display omics Item
    # for each omics data type
    # and for each dataser
    # if no error message
    observeEvent(inputData$loadData, {

      #continue only if message is true
      validate({
        need(validate.status == 0 ,message="ok")
      })

          #### Item for each omics #####
          output$omics <- renderMenu({
            menu_list <- list()
            menu_list <- list(
              menu_list,
              sidebarMenu(id = "sbm",

                  lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){

                    do.call(what = menuItem,
                        args = c(text = paste0(omics, " Analysis"), tabName = paste0(omics, "Analysis"),

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
                                  }
                                 )
                          )
                        )
                    })
                  )
              )
            sidebarMenu(.list = menu_list)
            })


    ##########################################
    # Part3 : Data Exploratory
    ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){

      switch(omics ,
             "RNAseq"={
                 lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){

                   callModule(RNAseqDataExplorTab, paste0("RNAseq",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
                  })
             },
             "proteomics"={
                 lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){

                   callModule(ProtMetaDataExplorTab, paste0("proteomics",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
                 })
             },
             "metabolomics"={
                 lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){

                   callModule(ProtMetaDataExplorTab, paste0("metabolomics",i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
                 })
              }
      )
    })

    ##########################################
    # Part4 : Data processing : filtering, Normalisation...
    ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){

      switch(omics ,
             "RNAseq"={

      lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
      #for(i in names(FlomicsMultiAssay@metadata$omicList[[omics]])){
        callModule(RNAseqDataNormTab, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
        #inputNorm[[paste0(omics, i)]] <- callModule(RNAseqDataNormTab, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
        
      #}
      })
             })
    })

    ##########################################
    # Part5 :  Diff Analysis
    ##########################################

    inputDiff <- list()
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){

      lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){
        
       # observeEvent(inputNorm[[paste0(omics, i)]]$validContrast, {
          
        inputDiff[[paste0(omics, i)]] <<- callModule(DiffExpAnalysis, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
      })
    })

  ##########################################
  # Part6 : Co-Expression Analysis
  ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){

      lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){

        observeEvent(inputDiff[[paste0(omics, i)]]$validContrast, {

          callModule(CoSeqAnalysis, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
        }, ignoreInit = TRUE)
      })
    })
    ##########################################
    # Part7 : Enrichment Analysis
    ##########################################
    lapply(names(FlomicsMultiAssay@metadata$omicList), function(omics){

      lapply(names(FlomicsMultiAssay@metadata$omicList[[omics]]), function(i){

        observeEvent(inputDiff[[paste0(omics, i)]]$validContrast, {

          callModule(AnnotationEnrichment, paste0(omics, i), FlomicsMultiAssay@metadata$omicList[[omics]][[i]])
        })
      })
    })


})


    ##########################################
    # Part8 : RMD REPORT
    ##########################################


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

      # Set up parameters to pass to Rmd document
      params <- list( FEdata = file.path(tempdir(), "FlomicsMultiAssay.RData"),
                      pngDir = tempdir())

      print(tempdir())
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))

      # rmarkdown::render(tempReport, output_file = file,
      #                   params = list( FEdata = file.path(tempdir(), "FlomicsMultiAssay.RData"),
      #                                  pngDir = tempdir()),
      #                   envir = new.env(parent = globalenv()))
    }
  )

    # # Automatically bookmark every time an input changes
    # observe({
    #   reactiveValuesToList(input)
    #   session$doBookmark()
    # })
    # # Update the query string
    # onBookmarked(updateQueryString)
})

