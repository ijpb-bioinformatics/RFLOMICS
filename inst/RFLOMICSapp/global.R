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
library(qgraph)
library(colourpicker)  ## ADDED August 2022 button to select colors
# library(circlize)
library(RColorBrewer) 
library(mixOmics) 
library(clusterProfiler) ## ADDED 27/10/22
library(org.At.tair.db) ## ADDED 27/10/22 # trouver un moyen de faire differemment. 
# library(KEGGREST) ## ADDED 27/10/22 # est-ce que ca marche ? Non, il faut kegg.db, qui  n'existe plus
source("00_module_common.R")
source("01_module_load_data.R")
source("02_module_set_stat_model.R")
source("03_module_data_explor.R")
source("04_module_diff_analysis.R")
source("05_module_coexp_analysis.R")
source("06_module_annot_enrichment.R")
source("06_module_annot_enrichment_clusterProf.R")
source("07_module_data_integration_MOFA.R")
source("08_module_data_integration_mixOmics.R")
source("coverPage.R")

#enableBookmarking(store = "server")
omics.dic <- list()
omics.dic[["RNAseq"]][["variableName"]] <- "genes"
omics.dic[["RNAseq"]][["valueType"]]    <- "counts"

omics.dic[["proteomics"]][["variableName"]] <- "proteins"
omics.dic[["proteomics"]][["valueType"]]    <- "XIC"

omics.dic[["metabolomics"]][["variableName"]] <- "metabolites"
omics.dic[["metabolomics"]][["valueType"]]    <- "XIC"


