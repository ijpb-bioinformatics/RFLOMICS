# library(testthat)
# library(RFLOMICS)

## Why is this code commented?
# Need either an internet connection or the file to test... 
# or org.at.tair

# ---- Construction MAE RFLOMICS ready for annotation analysis : ----
# data(ecoseed)
# # create rflomicsMAE object with ecoseed data
# MAE <- createRflomicsMAE(
#     projectName = "Tests",
#     omicsData   = list(ecoseed$RNAtest, ecoseed$metatest, ecoseed$protetest),
#     omicsNames  = c("RNAtest", "metatest", "protetest"),
#     omicsTypes  = c("RNAseq","metabolomics","proteomics"),
#     ExpDesign   = ecoseed$design,
#     factorRef   = ecoseed$factorRef)
# names(MAE) <- c("RNAtest", "metatest", "protetest")
# 
# formulae <- generateModelFormulae( MAE) 
# MAE <- setModelFormula(MAE, formulae[[1]])
# 
# contrastList <- generateExpressionContrast(object = MAE) |> 
#     purrr::reduce(rbind) |>
#     dplyr::filter(contrast %in% c("(temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS)",
#                                   "((temperatureLow_imbibitionEI - temperatureLow_imbibitionDS) + (temperatureMedium_imbibitionEI - temperatureMedium_imbibitionDS) + (temperatureElevated_imbibitionEI - temperatureElevated_imbibitionDS))/3",
#                                   "((temperatureElevated_imbibitionEI - temperatureLow_imbibitionEI) - (temperatureElevated_imbibitionDS - temperatureLow_imbibitionDS))" ))
# MAE <- MAE |>
#     setSelectedContrasts(contrastList) |>
#     filterLowAbundance(SE.name = "RNAtest")                          |>
#     runNormalization(SE.name = "RNAtest", normMethod = "TMM")        |>
#     runDiffAnalysis(SE.name = "RNAtest", method = "edgeRglmfit")     |>
#     runCoExpression(SE.name = "RNAtest",
#                     K = 2:12,
#                     replicates = 2) 

# ---- Annotation test function - DiffExpEnrichment ----
# 

#### Need org.At.tair.db to run
# test_that("it's running from diffExpAnal - GO - RNASeq", {
#   
#   # Selecting only one contrast
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, nameList = "(temperatureMedium - temperatureLow) in imbibitionDS" , SE.name = "RNAtest", database = "GO",
#                                    list_args = list(OrgDb = "org.At.tair.db", 
#                                                     keyType = "TAIR", 
#                                                     pvalueCutoff = 0.05),
#                                    domain = c("BP", "MF", "CC"))
#   })
#   
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "(temperatureMedium - temperatureLow) in imbibitionDS" , database = "GO", domain = "BP")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
#   
#   # All contrasts
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "GO",
#                                    list_args = list(OrgDb = "org.At.tair.db", 
#                                                     keyType = "TAIR", 
#                                                     pvalueCutoff = 0.05),
#                                    domain = c("BP", "MF", "CC"))
#     
#   })
#   
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS", database = "GO", domain = "BP")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
#   
#   
#   
# })

# test_that("it's running from diffExpAnal - Custom - RNASeq", {
#   
#   df_custom <- vroom::vroom(file = paste0(system.file(package = "RFLOMICS"), 
#                                           "/ExamplesFiles/ecoseed/AT_GOterm_EnsemblPlants.txt"),
#                             show_col_types = FALSE)
#   
#   MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "custom",
#                                  list_args = list(pvalueCutoff = 0.05),
#                                  col_term = "GO term accession", 
#                                  col_gene = "Gene stable ID",
#                                  col_name = "GO term name",
#                                  col_domain = "GO domain",
#                                  annot = df_custom)
#   
#   expect(!is.null(MAE[["RNAtest"]]@metadata$DiffExpEnrichAnal$custom$summary), failure_message = "Custom RNAseq didn't work (summary not present)")
#   
#   # Selecting only one contrast, custom file
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, nameList = "(temperatureMedium - temperatureLow) in imbibitionDS" , SE.name = "RNAtest", database = "custom",
#                                    list_args = list(pvalueCutoff = 0.05),
#                                    col_term = "GO term accession", 
#                                    col_gene = "Gene stable ID",
#                                    col_name = "GO term name",
#                                    col_domain = "GO domain",
#                                    annot = df_custom)
#   })
#   
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "(temperatureMedium - temperatureLow) in imbibitionDS" , database = "custom", domain = "biological_process")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")
#   
#   #### Need org.At.tair.db to run
#   # All contrasts, GO database
#   # expect_no_error({
#   #   MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "GO",
#   #                                  list_args = list(OrgDb = "org.At.tair.db", 
#   #                                                   keyType = "TAIR", 
#   #                                                   pvalueCutoff = 0.05),
#   #                                  domain = c("BP", "MF", "CC"))
#   #   
#   # })
#   
  # expect({
  #   obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS", database = "GO", domain = "BP")
  #   nrow(obj@result) > 0
  # }, failure_message = "(GO RNAseq from DiffExp) - There is no result in the enrichment metadata part.")

  # All contrasts, KEGG database
  # expect_no_error({
  #   MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", database = "KEGG",
  #                                  list_args = list(organism = "ath",
  #                                                   keyType = "kegg",
  #                                                   pvalueCutoff = 0.5))
  # 
  # })
  # 
  # expect({
  #   obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "(temperatureElevated - temperatureLow) in imbibitionDS", database = "KEGG")[["no-domain"]]
  #   nrow(obj@result) > 0
  # }, failure_message = "(KEGG RNAseq from DiffExp) - There is no result in the enrichment metadata part.")

# })


# ---- Annotation test function - CoExpression enrichment ----

#### Need org.At.tair.db to run
# test_that("it's running from CoExpAnal - GO - RNASeq", {
#   
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", from = "CoExpAnal", database = "GO",
#                                    list_args = list(OrgDb = "org.At.tair.db", 
#                                                     keyType = "TAIR", 
#                                                     pvalueCutoff = 0.05),
#                                    domain = c("BP", "MF", "CC"))
#     
#   })
#   
#   
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "cluster.1", database = "GO", domain = "BP", from = "CoExpAnal")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO RNAseq from CoExp) - There is no result in the enrichment metadata part.")
#   
#   
#   expect_no_error({
#     MAE <- runAnnotationEnrichment(MAE, SE.name = "RNAtest", nameList = c("cluster.1", "cluster.2") ,
#                                    from = "CoExpAnal", database = "GO",
#                                    list_args = list(OrgDb = "org.At.tair.db", 
#                                                     keyType = "TAIR", 
#                                                     pvalueCutoff = 0.05),
#                                    domain = c("BP", "MF", "CC"))
#     
#   })
#   
#   expect({
#     obj <- RFLOMICS:::getEnrichRes(MAE[["RNAtest"]], contrast = "cluster.1", database = "GO", domain = "BP", from = "CoExpAnal")
#     nrow(obj@result) > 0
#   }, failure_message = "(GO RNAseq from CoExp) - There is no result in the enrichment metadata part.")
#   
# })
