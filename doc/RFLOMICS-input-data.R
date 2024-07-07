## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  echo = FALSE,
  comment = "#>"
)


## ----echo=FALSE---------------------------------------------------------------
library(RFLOMICS)

## ----ImportDesign, output="asis"----------------------------------------------
data(ecoseed)
DT::datatable(ecoseed$design)

## ----RNAseq, output="asis"----------------------------------------------------
geneCount  <- ecoseed$RNAtest
DT::datatable(geneCount[1:10, 1:5])

## ----Proteomics, output="asis"------------------------------------------------
protAbundance <- ecoseed$protetest
DT::datatable(protAbundance[1:10, 1:5])

## ----Metabolomics, output="asis"----------------------------------------------
metaAbundance <- ecoseed$metatest
DT::datatable(metaAbundance[1:10, 1:5])

