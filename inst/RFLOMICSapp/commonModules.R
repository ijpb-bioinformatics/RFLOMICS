
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



