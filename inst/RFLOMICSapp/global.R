# Packages

library(magrittr)
library(edgeR)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(SummarizedExperiment)
library(FactoMineR)
library(ggplot2)
library(gridExtra)
library(MultiAssayExperiment)
library(reshape2)

step_tag <- 4
tmpDir <- "/Users/nbessoltane/rFlomics.report.tmp"
unlink(tmpDir, recursive=TRUE)

#if (file.exists(tmpDir)){
  
#} else {
  dir.create(tmpDir)
  #dir.create(file.path(tmpDir, "tmp"))
#}

