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

# Contrasts
Contrasts.List.Bidon <- data.frame(
            "idContrast" = paste("C",c(1:14),sep=""),
           "Contrasts" =
             c("WT.Treated - WT.Control",
             "dbMT.Treated - dbMT.Control",
             "oxMT.Treated - oxMT.Control",
             "siMT.Treated - siMT.Control",
              "dbMT.Control - WT.Control",
             "siMT.Control - WT.Control",
             "oxMT.Control - WT.Control",
             "dbMT.Control - siMT.Control",
             "dbMT.Control - siMT.Control",
             "dbMT.Treated - WT.Treated",
             "siMT.Treated - WT.Treated",
             "oxMT.Treated - WT.Treated",
             "dbMT.Treated - siMT.Treated",
             "dbMT.Treated - siMT.Treated"))

Contrasts.Coeff.bidon <- data.frame("C1"=c(0,0,0,1,0,0,0,0,0,0),
                                    "C2"=c(0,0,0,0,1,0,0,0,0,0))


