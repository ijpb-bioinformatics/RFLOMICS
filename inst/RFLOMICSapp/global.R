# Packages

library(plotly)
library(DT)
library(magrittr)
library(coseq)
library(edgeR)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(shinyWidgets)
library(SummarizedExperiment)
library(gridExtra)
library(MultiAssayExperiment)
library(FactoMineR)
library(ggplot2)
library(reshape2)
library(kableExtra)
library(shinyBS)
library(MOFA2) 
library(qgraph) ## ADDED August 2022
library(colourpicker)  ## ADDED August 2022 button to select colors
# library(circlize) ## ADDED August 2022
library(RColorBrewer) ## ADDED August 2022
library(mixOmics) ## ADDED 15/09/22
source("00_module_common.R")
source("01_module_load_data.R")
source("02_module_set_stat_model.R")
source("03_module_data_explor.R")
source("04_module_diff_analysis.R")
source("05_module_coexp_analysis.R")
source("06_module_annot_enrichment.R")
source("07_module_data_integration_MOFA.R")
source("08_module_data_integration_mixOmics.R") ### ADDED 15/09/2022
source("coverPage.R")

#enableBookmarking(store = "server")
omics.dic <- list()
omics.dic[["RNAseq"]][["variableName"]] <- "genes"
omics.dic[["RNAseq"]][["valueType"]]    <- "counts"

omics.dic[["proteomics"]][["variableName"]] <- "proteins"
omics.dic[["proteomics"]][["valueType"]]    <- "XIC"

omics.dic[["metabolomics"]][["variableName"]] <- "metabolites"
omics.dic[["metabolomics"]][["valueType"]]    <- "XIC"


