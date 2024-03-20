library(testthat)
library(RFLOMICS)

# ---- Construction MAE RFLOMICS ready for differential analysis : ----
MAE <- generateExample(
  annotation = FALSE,
  integration = FALSE
) 

# ---- Annotation test function - DiffExpEnrichment ----
# 

test_that("it's running from diffExpAnal - GO - RNASeq", {
  
  # Selecting only one contrast
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, nameList = "H1", SE.name = "RNAtest", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db", 
                                                    keyType = "TAIR", 
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H1", database = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  # All contrasts
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db", 
                                                    keyType = "TAIR", 
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))
    
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H2", database = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  
  
})

test_that("it's running from diffExpAnal - Custom - RNASeq", {
  
  df_custom <- vroom::vroom(file = paste0(system.file(package = "RFLOMICS"), 
                                          "/ExamplesFiles/GO_annotations/Arabidopsis_thaliana_Ensembl_55.txt"),
                            show_col_types = FALSE)
  
  MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "custom",
                                 list_args = list(pvalueCutoff = 0.05),
                                 col_term = "GO term accession", 
                                 col_gene = "Gene stable ID",
                                 col_name = "GO term name",
                                 col_domain = "GO domain",
                                 annot = df_custom)
  
  expect(!is.null(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$custom$summary), failure_message = "Custom RNAseq didn't work (summary not present)")
  
  # Selecting only one contrast, custom file
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, nameList = "H1", SE.name = "RNAtest", database = "custom",
                                   list_args = list(pvalueCutoff = 0.05),
                                   col_term = "GO term accession", 
                                   col_gene = "Gene stable ID",
                                   col_name = "GO term name",
                                   col_domain = "GO domain",
                                   annot = df_custom)
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H1", database = "custom", domain = "biological_process")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  
  # All contrasts, GO database
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db", 
                                                    keyType = "TAIR", 
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))
    
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H2", database = "GO", domain = "BP")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
  # All contrasts, KEGG database
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "KEGG",
                                   list_args = list(organism = "ath", 
                                                    keyType = "kegg", 
                                                    pvalueCutoff = 0.5))
    
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "H2", database = "KEGG")[["no-domain"]]
    nrow(obj@result) > 0
  }, failure_message = "(KEGG RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
  
})


# ---- Annotation test function - CoExpression enrichment ----

test_that("it's running from CoExpAnal - GO - RNASeq", {
  
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", from = "CoExpAnal", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db", 
                                                    keyType = "TAIR", 
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))
    
  })
  
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "cluster.1", database = "GO", domain = "BP", from = "CoExpAnal")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from CoExp) - There is no result in the enrichment metadata part.")
  
  
  expect_no_error({
    MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", nameList = c("cluster.1", "cluster.2") ,
                                   from = "CoExpAnal", database = "GO",
                                   list_args = list(OrgDb = "org.At.tair.db", 
                                                    keyType = "TAIR", 
                                                    pvalueCutoff = 0.05),
                                   domain = c("BP", "MF", "CC"))
    
  })
  
  expect({
    obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "cluster.1", database = "GO", domain = "BP", from = "CoExpAnal")
    nrow(obj@result) > 0
  }, failure_message = "(GO RNAseq from CoExp) - There is no result in the enrichment metadata part.")
  
})
