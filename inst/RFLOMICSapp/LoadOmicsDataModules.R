LoadOmicsDataUI <- function(id){


  ns <- NS(id)

  tagList(


    box(title = "Load omic data", status = "warning", width = 12, height = NULL,
        #h5("[warning] dataset with omics (type == none) was igored..."),
        fluidRow(
          column(2,
                 # omic type
                 selectInput(inputId=ns('omicType1'), label='Omics',
                             choices = c("None"="none", "RNAseq"="RNAseq",
                                         "Proteomics"="proteomics", "Metabolomics"="metabolomics"),
                             selected = "none")
          ),
          column(4,

                 # matrix count/abundance input
                 fileInput(inputId = ns("data1"), "omics count/abundance (Ex.)",
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          ),
          column(4,
                 # metadata/QC bioinfo
                 fileInput(inputId = ns("metadataQC1"), "QC or metadata (Ex.)",
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
          ),
          column(2,
                 # dataset Name
                 textInput(inputId=ns("DataName1"), label="Dataset name", value="set1")
          )

          #,
          #column(2, actionButton("removeFactor1", "",
          #                       icon=icon("times", class = NULL, lib = "font-awesome")))
        ),
        uiOutput(ns("toAddData2")),
        actionButton(inputId = ns("addData"),"Add data")
    ),
    actionButton(inputId = ns("loadData"),"load")
  )

}


LoadOmicsData <- function(input, output, session){


  # as soon as the "add data" buttom has been clicked
  # => a new select/file Input was display
  # =>
  addDataNum <- 1
  dataName.vec  <- c()
  observeEvent(input$addData, {

    # add input select for new data
    addDataNum <<- addDataNum + 1
    output[[paste("toAddData", addDataNum, sep="")]] <- renderUI({
      list(
        fluidRow(
            column(2,
                   # omic type
                   selectInput(inputId=session$ns(paste0('omicType', addDataNum)), label='Omics',
                               choices = c("None"="none", "RNAseq"="RNAseq",
                                           "Proteomics"="proteomics", "Metabolomics"="metabolomics"),
                               selected = "none")
                   ),
            column(4,
                   # matrix count/abundance input
                   fileInput(inputId=session$ns(paste0("data", addDataNum)), "Feature count/abundance (txt)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                   ),
            column(4,
                   # metadata/QC bioinfo
                   fileInput(inputId=session$ns(paste0("metadataQC", addDataNum)), "QC or metadata (txt)",
                             accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                   ),
            column(2,
                   # dataset Name
                   textInput(inputId=session$ns(paste0("DataName", addDataNum)), label="Dataset name", value=paste0("set", as.character(addDataNum)))
                   )
            ),
        uiOutput(session$ns(paste("toAddData",addDataNum + 1,sep="")))
        )
      })
    })

  ListOfmaeSlots <- list()
    # as soon as the "load data" buttom has been clicked
    # => selected input(s) are loaded with loadData() function
    # => create tmp dir for png report
    observeEvent(input$loadData, {

      #### load data
      ####
      print("# 5- Load omics data...")

      inputs <- list()

      SummarizedExperimentList <- list()
      omicList<- list()
      listmap <- list()

      #### list of omic data laoded from interface
      validate.status <<- 0
      dataName.vec <- c()
      for (k in 1:addDataNum){



        if(input[[paste0("omicType", k)]] != "none"){

          omicType <- input[[paste0("omicType", k)]]

          dataName <- paste0(omicType, ".", gsub("[[:space:]]", "", input[[paste0("DataName", k)]]))

          # check presence of dataname
          if(dataName == ""){
            showModal(modalDialog( title = "Error message", "Dataset names is required : dataset ", k ))
            validate.status <<- 1
          }

          # check duplicat dataset name
          if(any(duplicated(omicList[[omicType]])) == TRUE){
            showModal(modalDialog( title = "Error message", "Dataset names must be unique : dataset ", (1:addDataNum)[duplicated(dataName.vec)] ))
            validate.status <<- 1
          }


          inputs[[dataName]][["omicType"]] <- omicType  ####

          dataFile <- input[[paste0("data", k)]]
          inputs[[dataName]][["dataFile"]] <- input[[paste0("data", k)]]$datapath ####

          # check omics data
          if(is.null(dataFile)){
            showModal(modalDialog( title = "Error message", "omics counts/abundances matrix is required : dataset ", k ))
            validate.status <<- 1
          }

          qcFile   <- input[[paste0("metadataQC", k)]]
          inputs[[dataName]][["qcFile"]] <- input[[paste0("metadataQC", k)]]$datapath ####


          validate({
            need(validate.status == 0, message="error")
          })


        }
      }

      FlomicsMultiAssay <<- FlomicsMultiAssay.constructor(inputs = inputs, Design=Design)

    })

  return(input)

}

