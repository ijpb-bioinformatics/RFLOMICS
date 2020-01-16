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

