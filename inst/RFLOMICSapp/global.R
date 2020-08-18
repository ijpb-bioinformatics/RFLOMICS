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
source("DataExploratoryModules.R")
source("NormalizationModules.R")
source("DiffExpressionModules.R")
source("CoExpressionModules.R")
source("commonModules.R")



SupportedOmics <- c("RNAseq", "proteomics", "metabolomics")
