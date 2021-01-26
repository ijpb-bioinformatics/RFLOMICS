library(testthat)
library(RFLOMICS)

test_check("RFLOMICS")

test <- read.table("inst/ExamplesFiles/experimental_design.txt",header = TRUE,row.names = 1)
ExNames <- names(test)
ExTypes <- c("Bio","Bio","batch")

test1 <- GetModelFormulae(ExNames,ExTypes)
