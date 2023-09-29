
# shiny
library(shiny)
library(shinyWidgets)
library(shinyBS)

# data management
library(magrittr)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(reshape2)
library(kableExtra)
library(DT)

# graph / table
library(plotly)
library(gridExtra)
library(ggplot2)
# library(qgraph)
# library(colourpicker)  
library(circlize)
# library(RColorBrewer)
library(grid)
library(png)

# stat
library(edgeR)
library(FactoMineR)
library(coseq)

# parallel
library(parallel) # TODO change for biocparallel 

# #############
# MODULES

source("00_module_common.R")
source("01_module_load_data.R")
source("02_module_set_stat_model.R")
source("03_module_data_explor.R")
source("04_module_diff_analysis.R")
source("05_module_coexp_analysis.R")
source("06_module_annot_enrichment.R")
source("07_module_data_integration_MOFA.R")
source("08_module_data_integration_mixOmics.R") 
source("coverPage.R")


