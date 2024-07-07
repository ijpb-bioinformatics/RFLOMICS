## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(RFLOMICS)
#  runRFLOMICS()

## ----fig.show="hold", out.width="40%", fig.align = "center", echo=FALSE-------
knitr::include_graphics(c("Images/00_completeDesign.png","Images/00_nonBalancedDesign.png"))

## -----------------------------------------------------------------------------
# geneToKeep <- counts[rowSums(edgeR::cpm(counts) >= CPM_Cutoff) >= x, ]

## -----------------------------------------------------------------------------
# dge <- edgeR::DGEList(counts=counts, group=groups)
# dge <- edgeR::calcNormFactors(dge, method="TMM")

