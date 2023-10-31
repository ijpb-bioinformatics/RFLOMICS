library(RFLOMICS)
library(testthat)

# ---- Creating files to test ----

# # Create files with missing values, special char, etc.
# condEcoseed <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt"))
# set.seed(2000)
# modDelete1 <- sample(x = 1:nrow(condEcoseed), size = 3)
# set.seed(2001)
# modDelete2 <- sample(x = 1:nrow(condEcoseed), size = 3)
# # set.seed(2002)
# # modDelete3 <- sample(x = 1:nrow(condEcoseed), size = 3)
# 
# condEcoseedNA <- condEcoseed
# condEcoseedNA$Repeat[modDelete1] <- NA
# condEcoseedNA <- condEcoseedNA %>% tibble::rownames_to_column(var = "sample") %>% dplyr::relocate(sample)
# write.table(condEcoseedNA, file = "inst/testFiles/EcoseedCondNA.txt", sep = "\t", row.names = FALSE)
# 
# condEcoseedspecChar <- condEcoseed
# levels(condEcoseedspecChar$temperature) <- c("Ele_vated", "L&w", "Médium")
# condEcoseedspecChar <- condEcoseedspecChar %>% tibble::rownames_to_column(var = "sample") %>% dplyr::relocate(sample)
# write.table(condEcoseedspecChar, file = "inst/testFiles/EcoseedCondSpecChar.txt", sep = "\t", row.names = FALSE)

## NEED 
# TODO
# Create file with duplicated colnames
# Create file with duplicated modalities
# Same work with data

# ---- Tests read_exp_design ----

expect_no_error(read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt")))
expect_error(read_exp_design())

# Missing values
condEcoseedTest <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/testFiles/EcoseedCondNA.txt"))
# TODO verifier le constructeur ET lancer une analyse diff avec des NA

# Special characters
condEcoseedTest <- read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/testFiles/EcoseedCondSpecChar.txt"))
# TODO verifier le constructeur ET lancer une analyse diff avec des characteres speciaux

# ---- Tests read_omics_data ----
 
expect_no_error(RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")))
expect_error(RFLOMICS::read_omics_data(file = "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))

# ---- Tests RFLOMICS constructor ----

# load data

RNAdat <-  RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt"))


MAE <- RFLOMICS::FlomicsMultiAssay.constructor(
  list("RNAtest"     = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
                            "omicType" = "RNAseq"),
       "metatest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")), 
                         "omicType" = "metabolomics"),
       "protetest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")), 
                          "omicType" = "proteomics")
  ),
  projectName = "Tests", 
  ExpDesign = RFLOMICS::read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt")),
  refList = c("Repeat" = "rep1", "temperature" = "Low", "imbibition" = "DS"),
  typeList = c("Repeat" = "batch", "temperature" = "Bio", "imbibition" = "Bio"))


# TODO add test when several names are the same
MAE <- RFLOMICS::FlomicsMultiAssay.constructor(
  list("RNAtest"     = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/transcriptome_ecoseed.txt")),
                            "omicType" = "RNAseq"),
       "metatest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/metabolome_ecoseed.txt")),
                         "omicType" = "metabolomics"),
       "protetest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")),
                          "omicType" = "proteomics"),
       "protetest" = list("data" = RFLOMICS::read_omics_data(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/proteome_ecoseed.txt")),
                          "omicType" = "proteomics")
  ),
  projectName = "Tests",
  ExpDesign = RFLOMICS::read_exp_design(file = paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/condition.txt")),
  refList = c("Repeat" = "rep1", "temperature" = "Low", "imbibition" = "DS"),
  typeList = c("Repeat" = "batch", "temperature" = "Bio", "imbibition" = "Bio"))
