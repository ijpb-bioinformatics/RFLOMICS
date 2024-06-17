### --- Global import ----
#' @importFrom MOFA2 views_names get_factors plot_factor get_weights
#' plot_weights plot_variance_explained plot_factor_cor
#' plot_data_overview plot_data_heatmap get_dimensions

# ---- remove batch effects from omics : ----

#' @title .rbeFunction
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param SEobject An object of class \link{RflomicsSE}
#' @return An object of class \link{RflomicsSE}
#' @importFrom limma removeBatchEffect
#' @importFrom stats model.matrix as.formula
#' @keywords internal
#' @noRd
#'
.rbeFunction <- function(object, SEobject, cmd = FALSE) {
    assayTransform <- assay(SEobject)
    colBatch <- getBatchFactors(SEobject)
    
    if (cmd) {
        message(
            "#     =>Correction for Batch: ",
            paste(colBatch, collapse = " "),
            " in ",
            SEobject@metadata$omicType
        )
    }
    
    newFormula <-
        gsub(pattern = paste(paste(colBatch, "[+]"),  collapse = "|"),
             "",
             getModelFormula(SEobject))
    newFormula <- gsub(pattern = "~ [+] ", "~ ", newFormula) # use ?
    colData <- getDesignMat(SEobject)
    designToPreserve <-
        model.matrix(as.formula(newFormula), data = colData)
    
    if (length(colBatch) == 1) {
        rbeRes <-
            removeBatchEffect(assayTransform, batch = colData[, colBatch],
                              design = designToPreserve)
    } else if (length(colBatch) >= 2) {
        rbeRes <- removeBatchEffect(
            assayTransform,
            batch  = colData[, colBatch[1]],
            batch2 = colData[, colBatch[2]],
            design = designToPreserve
        )
    }
    if (length(colBatch) > 2)
        message("sorry, only 2 batches effect for now!")
    
    SEobject@metadata[["correction_batch_method"]] <-
        "limma (removeBatchEffect)"
    
    assay(SEobject) <- t(scale(t(rbeRes), center = TRUE, scale = FALSE))
    # assay(SEobject) <- rbeRes
    
    return(SEobject)
}

# ----- Transform rnaseq assay from SE ----

#' @title Remove batch effect and transform rnaseq data
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param SEname the name of the rnaseq data to transform.
#' Supposed to be a RflomicsSE.
#' @param correctBatch if TRUE, correction of batch effects.
#' @param transformation the name of the transformation to be applied
#' on counts. Default is limma voom.
#' @param variableLists vector of variable names
#' No other option for now.
#' @return An object of class \link{RflomicsMAE}
#' @importFrom stats model.matrix formula
#' @importFrom dplyr filter
#' @importFrom limma voom
#' @importFrom edgeR DGEList
#' @keywords internal
#' @noRd
#'

.rnaseqRBETransform <- function(object,
                                SEname,
                                correctBatch = FALSE,
                                transformation = "limma (voom)",
                                variableNames = NULL,
                                cmd = FALSE) {
    if (!is(object, "RflomicsMAE"))
        stop("object is not a RflomicsMAE")
    
    rnaDat <- object[[SEname]]
    assayTransform <- assay(rnaDat)
    
    if (!is.integer(assayTransform) &&
        !identical(assayTransform, floor(assayTransform))) {
        message("You indicated RNASeq data for ",
                SEname,
                "but it is not recognized as count data")
    }
    
    DMat      <- getDesignMat(rnaDat)
    coefNorm  <- getCoeffNorm(rnaDat)
    designMat <-
        model.matrix(formula(getModelFormula(rnaDat)), data = DMat)
    
    DGEObject <- DGEList(
        counts       = assayTransform,
        norm.factors = coefNorm$norm.factors,
        lib.size     = coefNorm$lib.size,
        samples      = DMat
    )
    
    limmaRes <- voom(DGEObject, design = designMat)
    
    assay(rnaDat) <- limmaRes$E
    
    if (correctBatch)
        rnaDat <- .rbeFunction(object, rnaDat)
    rnaDat <- rnaDat[variableNames,]
    
    rnaDat@metadata[["correction_batch"]]             <-
        correctBatch
    rnaDat@metadata[["transform_results_all"]]        <- limmaRes
    rnaDat@metadata[["transform_method_integration"]] <-
        transformation
    
    object[[SEname]] <- rnaDat
    
    return(object)
}


# ----- Transform rnaseq assay from SE ----

#' @title .rbeTransform
#'
#' @param object An object of class \link{RflomicsMAE}
#' @param SEname the name of the omics data to transform. No counts data.
#' @param correctBatch if TRUE, correction of batch effects.
#' @return An object of class \link{RflomicsMAE}
#' @keywords internal
#' @noRd
#'

.rbeTransform <- function(object,
                          SEname,
                          correctBatch = TRUE,
                          variableNames = NULL,
                          type = "union",
                          choice = "DE",
                          cmd = FALSE) {
    omicsDat <- object[[SEname]]
    metadata(omicsDat)[["correction_batch"]]             <-
        correctBatch
    metadata(omicsDat)[["transform_method_integration"]] <-
        getTransSettings(omicsDat)$method
    
    if (correctBatch)
        omicsDat <- .rbeFunction(object, omicsDat)
    omicsDat <- omicsDat[variableNames,]
    
    object[[SEname]] <- omicsDat
    
    return(object)
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

.plot_MO_1 <- function(Data_res) {
    dat_explained <- melt(do.call("rbind", Data_res$prop_expl_var))
    colnames(dat_explained) <- c("Dataset", "Component",
                                 "% of explained variance")
    dat_explained$`% of explained variance` <-
        dat_explained$`% of explained variance` * 100
    
    dat_comb <- dat_explained %>%
        group_by(Dataset) %>%
        summarise("Cumulative Explained Variance" = sum(`% of explained variance`))
    
    if (is(Data_res, "block.splsda") ||
        is(Data_res, "block.plsda")) {
        dat_comb <- dat_comb %>% filter(Dataset != "Y")
    }
    
    gg1 <- ggplot(dat_comb,
                  aes(x = Dataset,
                      y = `Cumulative Explained Variance`)) +
        geom_col(fill = "darkblue") +
        theme_bw() +
        theme(
            axis.text  = element_text(size = 12),
            axis.line  = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size = 12)
        ) +
        ylab("") + xlab("") +
        ggtitle("Cumulative explained variance")
    
    return(gg1)
}

#' @keywords internal
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_tile aes theme_classic theme element_blank
#' element_text scale_fill_gradientn ylab ggtitle
#' @importFrom dplyr filter
#' @noRd
.plot_MO_2 <- function(Data_res) {
    dat_explained <- melt(do.call("rbind", Data_res$prop_expl_var))
    colnames(dat_explained) <- c("Dataset", "Component",
                                 "% of explained variance")
    
    dat_explained$`% of explained variance` <-
        dat_explained$`% of explained variance` * 100
    
    if (is(Data_res, "block.splsda") ||
        is(Data_res, "block.plsda")) {
        dat_explained <- dat_explained %>% filter(Dataset != "Y")
    }
    
    # Chunk of code to be cohesive with MOFA2::plot_explained_variance
    gg2 <- ggplot(dat_explained, aes(x = Dataset, y = Component)) +
        geom_tile(aes(fill = `% of explained variance`), color = "gray88") +
        geom_text(aes(label = round(`% of explained variance`,2)), 
                  color = "white", size = 4) +
        theme_classic() +
        theme(
            axis.text  = element_text(size = 10),
            axis.line  = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size = 10),
        ) +  ylab("") + xlab("") +
        scale_fill_gradientn(
            name = "% of explained\nvariance",
            colors = c("gray97", "darkblue"),
            guide = "colorbar",
            limits = c(
                min(dat_explained$`% of explained variance`),
                max(dat_explained$`% of explained variance`)
            )
        ) +
        ggtitle("Percentage of explained variance \n per component per block")
    
    return(gg2)
}

######################## INTEGRATION USING MOFA ########################

#' @title Run MOFA Analysis
#' @description Runs a MOFA analysis based on an untrained MOFA object and
#' user arguments.
#' @param object An untrained MOFA object
#' @param scale_views boolean.
#' MOFA option to scale the views so they have the same variance.
#' Default is FALSE.
#' @param maxiter integer. MOFA option, maximum number of iterations
#' to be considered if there it does not converge. Default is 1000.
#' @param num_factors integer. MOFA option,
#' maximum number of factor to consider. Default is 10.
#' @param ... Not in use.
#' @return A list with an untrained MOFA object
#' (containing all options for the run) and a trained MOFA object
#' @importFrom MOFA2 get_default_data_options get_default_model_options
#' get_default_training_options prepare_mofa run_mofa
#' @keywords internal
#' @noRd
.runMOFAAnalysis <- function(object,
                             scale_views = FALSE,
                             maxiter = 1000,
                             num_factors = 10,
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

    MOFA.messages <- list()
    withCallingHandlers({
        MOFAObject.untrained <- prepare_mofa(
            object           = object,
            data_options     = data_opts,
            model_options    = model_opts,
            training_options = train_opts
        )
    }, warning = function(warn){
        MOFA.messages[[length(MOFA.messages) + 1]] <<- warn
    } 
    )
    outfile <- file.path(tempdir(), 
                         paste0("mofa_", 
                                format(Sys.time(), format = "%Y%m%d-%H%M%S"), 
                                ".hdf5"))
    MOFAObject.trained <-
        run_mofa(MOFAObject.untrained,
                 use_basilisk = FALSE,
                 outfile = outfile,
                 save_data = TRUE)
    
    return(
        list(
            "MOFAObject.untrained" = MOFAObject.untrained,
            "MOFAObject.trained"   = MOFAObject.trained,
            "MOFA.messages"        = MOFA.messages
        )
    )
}


######################## INTEGRATION USING MixOMICS ########################

#' @title Run MixOmics Analysis
#' @description Run a mixOmics analysis. Given the specification of the user
#' (type of response, if response there is, multi block or not)
#'  the function will determine which mixOmics function to use.
#'  Please see mixOmics manual or website for more information.
#' @param object list of blocks (matrices) with the same samples (rows)
#' @param scale_views Boolean. Do the matrices have to be scaled?
#' Is used inside the mixOmics function.
#' @param selectedResponse Default is NULL (pls functions).
#' The response, can be a matrix or a single factor
#' (discriminant analysis is set in this case).
#' @param ncomp Number of component to be computed.
#' @param link_dataset Numeric between 0 and 1.
#' Only used for multi block analysis.
#' Indicates the correlation to be considered between the matrices.
#' It impacts the weights of the features, hence the feature selection.
#' Please see mixOmics user's guide for better explanation.
#' @param link_response Numeric between 0 and 1.
#' Indicates the correlation to be considered between the matrices
#' and the response matrix.
#' It impacts the weights of the features, hence the feature selection.
#' Please see mixOmics user's guide for better explanation.
#' @param sparsity Boolean. Indicates if there is a feature selection purpose.
#' If TRUE, functions like spls(da), block.spls(da) will be used.
#' @param cases_to_try If sparsity is TRUE,
#' indicates the number of cases to try for the feature selection.
#'  The best outcome, as computed by tuning function, will be displayed.
#' @return A mixOmics result.
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr arrange select select_if
#' @importFrom tidyselect all_of
#' @importFrom mixOmics plsda splsda block.plsda block.splsda tune.block.splsda
#' tune.splsda unmap
#' @importFrom utils getFromNamespace
#' @keywords internal
#' @noRd

.runMixOmicsAnalysis <- function(object,
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
    
    if (length(intersect(colnames(object$metadata), selectedResponse)) == 0) {
        stop("Selected Responses are not columns from metadata")
    }
    
    Y <- data.frame(object$metadata, stringsAsFactors = TRUE)
    Y <- Y[, selectedResponse, drop = FALSE]
    
    # Is this a discriminant analysis?
    # TODO : is this used?
    if (ncol(Y) == 1 && is.factor(Y[, 1])) {
        dis_anal <- TRUE
        Y <- Y[, 1]
    } else {
        YrowNames <- rownames(Y)
        YFactors <-
            do.call("cbind", lapply(
                seq_len(ncol(Y)),
                FUN = function(j) {
                    if (is.factor(Y[, j])) {
                        mat_return <- unmap(Y[, j])
                        colnames(mat_return) <-
                            attr(mat_return, "levels")
                        return(mat_return)
                    }
                }
            ))
        
        Y <- cbind(Y %>% select_if(is.numeric), YFactors)
        Y <- apply(Y, 2, as.numeric)
        rownames(Y) <- YrowNames
    }
    
    # Design matrix
    Design_mat <- matrix(
        link_datasets,
        nrow = length(object$blocks) + 1,
        ncol = length(object$blocks) + 1
    )
    Design_mat[, ncol(Design_mat)] <-
        Design_mat[nrow(Design_mat), ] <- link_response
    diag(Design_mat) <- 0
    
    # What function to use for the analysis
    # (can't be sparse if there is a continous response)
    functionName <- "pls"
    if (dis_anal)
        functionName <- paste0(functionName, "da")
    if (sparsity &&
        !is.numeric(Y))
        functionName <- paste0("s", functionName)
    if (length(object$blocks) > 1)
        functionName <- paste0("block.", functionName)
    
    # Model Tuning (if required, for sparsity)
    if (sparsity && dis_anal && functionName != "block.spls") {
        tune_function <- paste0("tune.", functionName)
        
        test_keepX <- lapply(
            object$blocks,
            FUN = function(dat) {
                ceiling(seq(
                    from = ceiling(0.1 * ncol(dat)),
                    to = ncol(dat),
                    length.out = cases_to_try
                ))
            }
        )
        
        ResponseY <<- Y
        
        list_tuning_args <- list(
            X = object$blocks,
            Y = Y,
            ncomp = ncomp,
            scale = scale_views,
            test.keepX = test_keepX,
            folds = min(table(Y) - 1, 10)  
        )
        if (length(object$blocks) > 1)
            list_tuning_args$design <- Design_mat
        
        list_res$tuning_res <-
            do.call(getFromNamespace(tune_function, ns = "mixOmics"),
                    list_tuning_args)
    }
    
    # Model fitting
    list_args <- list(Y = Y,
                      ncomp = ncomp,
                      scale = scale_views)
    
    if (length(object$blocks) == 1) {
        list_args$X <- object$blocks[[1]]
    } else {
        list_args$X <- object$blocks
    }
    if (length(object$blocks) > 1)
        list_args$design <- Design_mat
    if (sparsity &&
        !functionName %in% c("block.spls", "block.pls")) {
        list_args$keepX <- list_res$tuning_res$choice.keepX
    }
    
    list_res$analysis_res <-
        do.call(getFromNamespace(functionName, ns = "mixOmics"),
                list_args)
    
    return(list_res)
}


# ---- MOFA relationship to factor ----

#' @title Run test to compute relationship between MOFA2 factors and features
#' @description 
#' Run a kruskal.test between all factors and all biological and metadata 
#' factors entered by the user in the interface. Is used inside the interface
#' and the report. 
#' @param mofaRes an object from a MOFA run. 
#' @param method P value adjustment method. One of "BH", "Bon", etc. Same as 
#' in p.adjust method. 
#' @param ... Not in use at the moment. 
#' @return A table or a graph. 
#' @importFrom stats kruskal.test
#' @importFrom MOFA2 get_factors
#' @keywords internal
#' @noRd

.relationsMOFA <- function(mofaRes, method = "BH", ...){
    factors   <- get_factors(mofaRes)
    ExpDesign <- mofaRes@samples_metadata
    ExpDesign$group  <- NULL
    ExpDesign$sample <- NULL
    
    res_aov <- lapply(
        seq_len(ncol(ExpDesign)),
        FUN = function(i) {
            p.adjust(unlist(lapply(seq_len(ncol(factors$group1)),
                          FUN = function(j){
                              kruskal.test(x = factors$group1[,j], 
                                           g = ExpDesign[,i])$p.value
                          })), method = method)
        }
    )
    names(res_aov) <- colnames(ExpDesign)
    
    res_res <- data.frame(do.call("rbind", res_aov))
    # colnames(res_res) <- gsub("Response ", "", colnames(res_res))
    colnames(res_res) <- paste0("Factor ", seq_len(ncol(res_res)))
    
    return(res_res)
}
