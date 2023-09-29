### --- Global import ----
#' @importFrom MOFA2 views_names get_factors plot_factor get_weights 
#' plot_weights plot_variance_explained plot_factor_cor
#' plot_data_overview plot_data_heatmap get_dimensions
 
# ---- remove batch effects from omics : ----

#' @title rbe_function
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEobject An object of class \link{SummarizedExperiment}
#' @return An object of class \link{SummarizedExperiment}
#' @importFrom limma removeBatchEffect
#' @importFrom stats model.matrix as.formula
#' @importFrom SummarizedExperiment assay
#' @export
#' @noRd
#'
rbe_function = function(object, SEobject, cmd = FALSE){
  
  assayTransform <- assay(SEobject)
  ftypes <- getFactorTypes(object)
  colBatch <- names(ftypes)[ftypes == "batch"]
  
  if (cmd) {
    print(paste0("#     =>Correction for Batch: ", 
                 paste(colBatch, collapse = " "), 
                 " in ", SEobject@metadata$omicType))
  }
  
  newFormula <- gsub(pattern = paste(paste(colBatch, "[+]"),  collapse = "|"), "", getModelFormula(object))
  newFormula <- gsub(pattern = "~ [+] ", "~ ", newFormula) # use ?
  designToPreserve <- model.matrix(as.formula(newFormula), data = SEobject@metadata$Groups)
  
  if (length(colBatch) == 1) {
    rbeRes <- removeBatchEffect(assayTransform, batch = SEobject@metadata$Groups[,colBatch], design = designToPreserve)
  }else if (length(colBatch) >= 2) {
    
    rbeRes <- removeBatchEffect(assayTransform,
                                batch  = SEobject@metadata$Groups[,colBatch[1]],
                                batch2 = SEobject@metadata$Groups[,colBatch[2]],
                                design = designToPreserve)
  }
  if (length(colBatch) > 2) print("sorry, only 2 batches effect for now!") 
  
  SEobject@metadata[["correction_batch_method"]] <- "limma (removeBatchEffect)"
  
  assay(SEobject) <- rbeRes
  
  return(SEobject)
}

# ----- Transform rnaseq assay from SE ----

#' @title Remove batch effect and transform rnaseq data
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEname the name of the rnaseq data to transform. Supposed to be a SummarizedExperiment.
#' @param correctBatch if TRUE, correction of batch effects. 
#' @param transformation the name of the transformation to be applied on counts. Default is limma voom. 
#' No other option for now. 
#' @return An object of class \link{MultiAssayExperiment}
#' @importFrom stats model.matrix formula
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr filter
#' @importFrom limma voom
#' @importFrom edgeR DGEList
#' @export
#' @noRd
#'

rnaseqRBETransform <- function(object, 
                               SEname, 
                               correctBatch = FALSE, 
                               transformation = "limma (voom)",
                               contrasts_names = NULL,
                               type = "union",
                               choice = "DE",
                               cmd = FALSE
){
  
  # TODO : check if differential analysis has been performed.
  # TODO : in case of removeBatchEffect not working, what do we do?
  # TODO : if no DE (ex: intersection is 0 genes), what do we do?
  
  if (!is(object, "MultiAssayExperiment")) stop("object is not a MultiAssyExperiment")
  
  rnaDat <- object[[SEname]] 
  assayTransform <- assay(rnaDat)
  
  if (!is.integer(assayTransform) && !identical(assayTransform, floor(assayTransform))) {
    message("You indicated RNASeq data for ", SEname, "but it is not recognized as count data") 
    print(assayTransform[1:3, 1:3])
  }
  
  DMat      <- getDesignMat(object)
  coefNorm  <- getNormCoeff(rnaDat)
  designMat <- model.matrix(formula(getModelFormula(object)), data = DMat)
  
  DGEObject <- DGEList(
    counts       = assayTransform,
    norm.factors = coefNorm$norm.factors,
    lib.size     = coefNorm$lib.size,
    samples      = DMat %>% 
      filter(row.names(DMat) %in% colnames(assayTransform))
  )
  
  limmaRes <- voom(DGEObject,
                   design = designMat[which(rownames(designMat) %in% colnames(assayTransform)),]) 
  
  assay(rnaDat) <- limmaRes$E
  
  if (correctBatch)    rnaDat <- rbe_function(object, rnaDat)
  
  if (choice == "DE")  rnaDat <- rnaDat[opDEList(rnaDat, contrasts = contrasts_names, operation = type),]
  
  rnaDat@metadata[["correction_batch"]]             <- correctBatch
  rnaDat@metadata[["transform_results_all"]]        <- limmaRes 
  rnaDat@metadata[["transform_method_integration"]] <- transformation
  
  object[[SEname]] <- rnaDat
  
  return(object)
} 


# ----- Transform rnaseq assay from SE ----

#' @title RBETransform
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param SEname the name of the omics data to transform. No counts data.  
#' @param correctBatch if TRUE, correction of batch effects. 
#' @return An object of class \link{MultiAssayExperiment}
#' @export
#' @noRd
#'

RBETransform <- function(object,
                         SEname,
                         correctBatch = TRUE,
                         contrasts_names = NULL,
                         type = "union",
                         choice = "DE",
                         cmd = FALSE){
  
  omicsDat <- object[[SEname]]
  omicsDat@metadata[["correction_batch"]]             <- correctBatch
  omicsDat@metadata[["transform_method_integration"]] <- omicsDat@metadata$transform_method
  
  if (correctBatch)    omicsDat <- rbe_function(object, omicsDat)
  if (choice == "DE")  omicsDat <- omicsDat[opDEList(omicsDat, contrasts = contrasts_names, operation = type),]
  
  object[[SEname]] <- omicsDat
  
  return(object)
}

# ----- MixOmics: plot variance explained ----

#' @title plot_MO_varExp
#'
#' @param object An object of class \link{MultiAssayExperiment}
#' @param selectedResponse a character string of the response variable to consider
#' @param mode Can be NULL (default), "cumulative" or "comp". Defines the type of graph to return
#' @return An object of class \link{MultiAssayExperiment}
#' @importFrom ggpubr ggarrange
#' @export
#'

plot_MO_varExp <- function(object, selectedResponse, mode = NULL){
  
  if (class(object) != "MultiAssayExperiment") stop("Object is not a MultiAssayExperiment.")
  if (is.null(object@metadata$mixOmics)) stop("It seems this object has no mixOmics results.")
  if (is.null(object@metadata$mixOmics[[selectedResponse]])) stop("It seems you didn't run MixOmics on this particular variable.")
  
  Data_res <- object@metadata$mixOmics[[selectedResponse]]$MixOmics_results
  gg_return <- NULL
  
  if (is.null(mode)) {
    gg_return <- ggarrange(plot_MO_1(Data_res),
                           plot_MO_2(Data_res), 
                           ncol = 2)
  }
  else if (tolower(mode) == "cumulative") {
    gg_return <- plot_MO_1(Data_res)
  }
  else if (tolower(mode) == "comp") {
    gg_return <- plot_MO_2(Data_res)
  }
  
  return(gg_return)
}

#' @title Plot cumulative Explained variance for mixOmics results
#'
#' @param object a MAE object (produced by Flomics).
#' @param selectedResponse a character string. 
#' @return A ggplot2 graph
#' @keywords internal
#' @importFrom dplyr group_by summarise filter
#' @importFrom reshape2 melt
#' @noRd

plot_MO_1 <- function(Data_res){
  dat_explained <- melt(do.call("rbind", Data_res$prop_expl_var))
  colnames(dat_explained) <- c("Dataset", "Component", "% of explained variance")
  dat_explained$`% of explained variance` <- dat_explained$`% of explained variance`*100
  
  dat_comb <- dat_explained %>% 
    group_by(Dataset) %>% 
    summarise("Cumulative Explained Variance" = sum(`% of explained variance`))
  
  if (is(Data_res, "block.splsda") || is(Data_res, "block.plsda")) {
    dat_comb <- dat_comb %>% filter(Dataset != "Y")
  }
  
  gg1 <- ggplot(dat_comb, aes(x = Dataset, y = `Cumulative Explained Variance`)) +
    geom_col(fill = "darkblue") + 
    theme_classic() +
    theme(
      axis.text  = element_text(size = 12),
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 12)) + 
    ylab("") + 
    ggtitle("Cumulative explained variance")  
  
  return(gg1)
}

#' @keywords internal
#' @importFrom reshape2 melt
#' @importFrom dplyr filter
#' @noRd
plot_MO_2 <- function(Data_res){
  dat_explained <- melt(do.call("rbind", Data_res$prop_expl_var))
  colnames(dat_explained) <- c("Dataset", "Component", "% of explained variance")
  dat_explained$`% of explained variance` <- dat_explained$`% of explained variance`*100
  
  if (is(Data_res, "block.splsda") || is(Data_res, "block.plsda")) {
    dat_explained <- dat_explained %>% filter(Dataset != "Y")
  }
  
  # Chunk of code to be cohesive with MOFA2::plot_explained_variance
  gg2 <- ggplot(dat_explained, aes(x = Dataset, y = Component)) +
    geom_tile(aes(fill = `% of explained variance`)) + 
    theme_classic() +
    theme(
      axis.text  = element_text(size = 12),
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 12),
    ) +  ylab("") +  
    scale_fill_gradientn(colors = c("gray97","darkblue"), guide = "colorbar", 
                         limits = c(min(dat_explained$`% of explained variance`),
                                    max(dat_explained$`% of explained variance`))) + 
    ggtitle("Percentage of explained variance \n per component per block")
  
  return(gg2)
}


#---- MOFA2 - plot network ----

#' @title Plot correlation network between omics from MOFA2 results
#'
#' @param object a MOFA result
#' @param factor_choice chose weight from this factor
#' @param abs_weight_network threshold of weight to select entities to print
#' @param abs_min_cor_network correlation threshold 
#' @param network_layout one of spring, circle , circle + omics
#' @param omics_colors named list of colors palettes, one for each omics (can be NULL) 
#' @param posCol colors of positive edges
#' @param negCol colors of negative edges
#' @importFrom MOFA2 get_factors get_weights
#' @importFrom stats cor
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate arrange desc group_by filter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 melt
#' @importFrom qgraph qgraph
#' @importFrom ggpubr ggarrange
#' @return A ggplot2 graph
#' @export
#' 
MOFA_cor_network <- function(resMOFA,
                             factor_choice = 1,
                             abs_weight_network = 0.5,
                             abs_min_cor_network = 0.5,
                             network_layout = "Circle + omics",
                             omics_colors = NULL,
                             posCol = "red",
                             negCol = "blue"
){
  
  if (!is(resMOFA, "MOFA")) {
    stop("resMOFA has to be a MOFA results")
  }
  
  # Correlation matrix is done on all ZW, not on the selected factor. 
  data_reconst_list <- lapply(get_weights(resMOFA), FUN = function(mat){
    get_factors(resMOFA)$group1 %*% t(mat)})
  data_reconst <- do.call(cbind, data_reconst_list)
  cor_mat <- cor(data_reconst)
  
  
  features_metadata <- do.call(rbind, lapply(1:length(get_weights(resMOFA)), FUN = function(i){
    mat_weights <- data.frame(get_weights(resMOFA, scale = TRUE)[[i]])
    mat_weights$Table <- names(get_weights(resMOFA))[i]
    return(mat_weights)
  }))
  
  factor_selected <- paste0("Factor", factor_choice)
  
  feature_filtered <- features_metadata %>% 
    rownames_to_column("EntityName") %>%
    mutate(F_selected = abs(get(factor_selected))) %>% 
    arrange(desc(abs(F_selected))) %>% 
    group_by(Table) %>% 
    filter(abs(F_selected) > abs_weight_network)
  
  if (nrow(feature_filtered) > 0) {
    
    if(is.null(omics_colors)){
      omics_colors <- as.list(rownames(brewer.pal.info)[1:length(unique(feature_filtered$Table))]) 
      names(omics_colors) <- unique(feature_filtered$Table)
    }
    
    omics_colors <- lapply(unique(feature_filtered$Table), FUN = function(omicTable){
      brewer.pal(name = omics_colors[[omicTable]], n = 5)
    })
    names(omics_colors) <- unique(feature_filtered$Table)
    
    feature_filtered <- feature_filtered %>% group_by(Table) %>%
      mutate(Color = cut(F_selected, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))
    feature_filtered$Color2 <- sapply(1:nrow(feature_filtered), 
                                      FUN = function(i) omics_colors[[feature_filtered$Table[i]]][as.numeric(feature_filtered$Color[i])])
    
    # Layout
    layout_arg <- tolower(network_layout)
    if (tolower(layout_arg) == tolower("Circle + omics")) {
      layout_arg <- "groups"
    }
    
    # Network main graph
    cor_display <- cor_mat[rownames(cor_mat) %in% feature_filtered$EntityName, colnames(cor_mat) %in% feature_filtered$EntityName]
    
    if (any(abs(cor_display[upper.tri(cor_display)]) >= abs_min_cor_network)) {
      qgraph_plot <- qgraph(cor_display, minimum = abs_min_cor_network, 
                            cut = 0,
                            shape = "rectangle", labels = rownames(cor_display), vsize2 = 2, 
                            vsize = sapply(rownames(cor_display), nchar)*1.1,  layout = layout_arg,
                            esize = 2,
                            groups = gsub("[.]filtred", "", features_metadata$Table[match(rownames(cor_display), rownames(features_metadata))]),
                            posCol = posCol, negCol = negCol, 
                            details = FALSE,  legend = FALSE,
                            color = feature_filtered$Color2[match(rownames(cor_display), feature_filtered$EntityName)])
      qgraph_plot <- recordPlot()
      
      # Legend
      legend_matrix <- do.call("rbind", omics_colors)
      colnames(legend_matrix) <- c("(0,0.2]", "(0.2,0.4]", "(0.4,0.6]", "(0.6,0.8]", "(0.8, 1]")
      rownames(legend_matrix) <- gsub("[.]filtred", "", rownames(legend_matrix))
      legend.reshape <- melt(legend_matrix)
      
      gg.legend <-  ggplot(legend.reshape, aes(x = Var2, y = Var1)) + 
        geom_tile(fill = legend.reshape$value) + xlab("") + ylab("") + 
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
                           axis.ticks.y = element_blank(),
                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      
      gg.legend <- ggarrange(gg.legend, nrow = 3, ncol = 1) 
      
      # Actual plotting
      ggarrange(plotlist = list(qgraph_plot, gg.legend), nrow = 1, ncol = 2, widths = c(3, 1))
      
      
      
    }else{
      print("There is nothing to plot. Please try to lower the absolute weight, the correlation threshold or change the factor.")
    }
  }else{
    print("There is nothing to plot. Please try to lower the absolute weight, the correlation threshold or change the factor.")
  }
}

######################## INTEGRATION USING MOFA ########################

#' @title Run MOFA Analysis
#' @description Runs a MOFA analysis based on an untrained MOFA object and user arguments.
#' @param object An untrained MOFA object
#' @param scale_views boolean. MOFA option to scale the views so they have the same variance. Default is FALSE.
#' @param maxiter integer. MOFA option, maximum number of iterations to be considered if there it does not converge. Default is 1000.
#' @param num_factors integer. MOFA option, maximum number of factor to consider. Default is 10.
#' @param ... Not in use.
#' @return A list with an untrained MOFA object (containing all options for the run) and a trained MOFA object
#' @importFrom MOFA2 get_default_data_options get_default_model_options get_default_training_options prepare_mofa run_mofa
#' @importFrom reticulate py_capture_output
#' @keywords internal
#' 
runMOFAAnalysis <- function(object,
                            scale_views = FALSE,
                            maxiter = 1000,
                            num_factors = 10,
                            silent = TRUE,
                            ...) {
  
  if (!is(object, "MOFA")) {
    stop("object has to be a MOFA results")
  }
  
  data_opts  <- get_default_data_options(object)
  model_opts <- get_default_model_options(object)
  train_opts <- get_default_training_options(object)
  
  data_opts$scale_views  <- scale_views
  train_opts$maxiter     <- maxiter
  train_opts$verbose     <- FALSE
  model_opts$num_factors <- num_factors
  
  if (silent) {
    MOFAObject.untrained <- suppressMessages(suppressWarnings(
      prepare_mofa(
        object           = object,
        data_options     = data_opts,
        model_options    = model_opts,
        training_options = train_opts
      )))
  } else {
    MOFAObject.untrained <- prepare_mofa(
      object           = object,
      data_options     = data_opts,
      model_options    = model_opts,
      training_options = train_opts
    )
  }
  
  if (silent && require(reticulate)) {
    
    pycapt <- py_capture_output(
      {
        MOFAObject.trained <- suppressWarnings(suppressMessages( 
          run_mofa(MOFAObject.untrained, use_basilisk = FALSE, save_data = TRUE)
        ))
      }
    )
  } else {
    MOFAObject.trained <- run_mofa(MOFAObject.untrained, use_basilisk = FALSE, 
                                   save_data = TRUE)
  }
  
  # peut poser probleme au niveau python et mofapy.
  # Installer python, numpy et mofapy, ensuite reinstaller totalement package MOFA2 et restart R.
  
  return(list("MOFAObject.untrained" = MOFAObject.untrained, 
              "MOFAObject.trained"   = MOFAObject.trained))
}


######################## INTEGRATION USING MixOMICS ########################

#' @title Run MixOmics Analysis
#' @description Run a mixOmics analysis. Given the specification of the user (type of response, if response there is, multi block or not)
#'  the function will determine which mixOmics function to use. Please see mixOmics manual or website for more information.
#' @param object list of blocks (matrices) with the same samples (rows)
#' @param scale_views Boolean. Do the matrices have to be scaled? Is used inside the mixOmics function.
#' @param selectedResponse Default is NULL (pls functions). The response, can be a matrix or a single factor (discriminant analysis is set in this case).
#' @param ncomp Number of component to be computed.
#' @param link_dataset Numeric between 0 and 1. Only used for multi block analysis. Indicates the correlation to be considered between the matrices.
#'        It impacts the weights of the features, hence the feature selection. Please see mixOmics user's guide for better explanation.
#' @param link_response Numeric between 0 and 1. Indicates the correlation to be considered between the matrices and the response matrix.
#'        It impacts the weights of the features, hence the feature selection. Please see mixOmics user's guide for better explanation.
#' @param sparsity Boolean. Indicates if there is a feature selection purpose. If TRUE, functions like spls(da), block.spls(da) will be used.
#' @param cases_to_try If sparsity is TRUE, indicates the number of cases to try for the feature selection. The best outcome, as computed by tuning function, will be displayed.
#' @return A mixOmics result. 
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr arrange select select_if
#' @importFrom tidyselect all_of
#' @importFrom mixOmics plsda splsda block.plsda block.splsda tune.block.splsda tune.splsda
#' unmap
#' @keywords internal 

runMixOmicsAnalysis <- function(object,
                                scale_views = FALSE,
                                selectedResponse = NULL,
                                ncomp = 2,
                                link_datasets = 1,
                                link_response = 1,
                                sparsity = FALSE,
                                cases_to_try = 5,
                                ...) {
  
  list_res <- list()
  dis_anal <- FALSE # is this a discriminant analysis
  
  # TODO : check this, c'est un peu etrange d'avoir eu a reordonner les lignes dans prepareForIntegration
  # du coup ca demande de le faire aussi pour Y
  # TODO : faudrait le faire directement dans preparedList ?
  Y <- data.frame(object$metadata, stringsAsFactors = TRUE) %>%
    rownames_to_column(var = "rowNam") %>%
    arrange(rowNam) %>%
    column_to_rownames(var = "rowNam") %>%
    select(all_of(selectedResponse))
  
  # Is this a discriminant analysis?
  # TODO : is this used?
  if (ncol(Y) == 1 && is.factor(Y[, 1])) {
    dis_anal <- TRUE
    Y <- Y[, 1]
  } else {
    YrowNames <- rownames(Y)
    YFactors <- do.call("cbind", lapply(1:ncol(Y), FUN = function(j) {
      if (is.factor(Y[, j])) {
        mat_return <- unmap(Y[, j])
        colnames(mat_return) <- attr(mat_return, "levels")
        return(mat_return)
      }
    }))
    
    Y <- cbind(Y %>% select_if(is.numeric), YFactors)
    Y <- apply(Y, 2, as.numeric)
    rownames(Y) <- YrowNames
  }
  
  # Design matrix
  Design_mat <- matrix(link_datasets,
                       nrow = length(object$blocks) + 1,
                       ncol = length(object$blocks) + 1
  )
  Design_mat[, ncol(Design_mat)] <-
    Design_mat[nrow(Design_mat), ] <- link_response
  diag(Design_mat) <- 0
  
  # What function to use for the analysis (can't be sparse if there is a continous response)
  functionName <- "pls"
  if (dis_anal) functionName <- paste0(functionName, "da")
  if (sparsity && !is.numeric(Y)) functionName <- paste0("s", functionName)
  if (length(object$blocks) > 1) functionName <- paste0("block.", functionName)
  
  # Model Tuning (if required, for sparsity)
  if (sparsity && dis_anal && functionName != "block.spls") {
    # no tune.block.spls so far, there is one for spls
    # It is a bit weird to consider tuning with folds when there is very few samples per condition
    # Add a warning or something when it's the case?
    # Also consider adding a warning when the number of feature is still very high even after differential analysis
    
    tune_function <- paste0("tune.", functionName)
    
    test_keepX <- lapply(object$blocks, FUN = function(dat) {
      ceiling(seq(from = ceiling(0.1 * ncol(dat)), to = ncol(dat), length.out = cases_to_try))
    })
    
    list_tuning_args <- list(
      X = object$blocks,
      Y = Y,
      ncomp = ncomp,
      scale = scale_views,
      test.keepX = test_keepX,
      folds = min(floor(nrow(Y) / 2), 10) # TODO ???
    )
    if (length(object$blocks) > 1) list_tuning_args$design <- Design_mat
    
    list_res$tuning_res <- do.call(getFromNamespace(tune_function, ns = "mixOmics"), list_tuning_args)
  }
  
  # Model fitting
  
  list_args <- list(
    X = object$blocks,
    Y = Y,
    ncomp = ncomp,
    scale = scale_views
  )
  if (length(object$blocks) > 1) list_args$design <- Design_mat
  if (sparsity && !functionName %in% c("block.spls", "block.pls")) list_args$keepX <- list_res$tuning_res$choice.keepX
  
  list_res$analysis_res <- do.call(getFromNamespace(functionName, ns = "mixOmics"), list_args)
  
  return(list_res)
}


