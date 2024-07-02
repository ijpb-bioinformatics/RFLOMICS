
################################### CO-EXPRESSION #############################

#' @title runCoExpression
#' @description This method performs a co-expression/co-abundance analysis of 
#' omic-data.
#' @details For now, only the coseq function of the coseq package is used.
#' For RNAseq data, parameters used are those recommended in DiCoExpress 
#' workflow (see the reference). This parameters are: \code{model="normal"}, 
#' \code{transformation="arcsin"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="TMM"}, \code{meanFilterCutoff = 50}
#' For proteomic or metabolomic, data are scaled by protein or metabolite 
#' to group them by expression profiles rather than by expression intensity.
#' After data scaling, recommended parameters (from \code{coseq} developers) 
#' for co-expression analysis are:
#' \code{model="normal"}, \code{transformation="none"}, 
#' \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="none"},  \code{meanFilterCutoff = NULL}.
#' @return
#' An S4 object of class \link{RflomicsSE}
#' All the results are stored as a named list \code{CoExpAnal}
#'  in the metadata slot of a given \code{RflomicsSE} object. Objects are:
#' The runCoExpression method return several results, for \link{coseq} 
#' method, objects are:
#' \itemize{
#' \item{\code{setting:} }{co-expression analysis settings. See \code{getCoexpSetting}}
#' \item{\code{results:} }{boolean indicating if the co-expression analysis succeed}
#' \item{\code{coseqResults:} }{the raw results of \code{coseq}}
#' \item{\code{clusters:} }{a List of clusters}
#' \item{\code{cluster.nb:} }{The number of cluster}
#' \item{\code{plots:} }{The plots of \code{coseq} results}
#' \item{\code{stats:} }{A tibble summarising failed jobs: reason, propoif any}
#' }
#' @param object An object of class \link{RflomicsSE}
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @param contrastNames names of the contrasts from which the DE entities 
#' have to be taken. Can be NULL, in that case every contrasts from the differential 
#' analysis are taken into consideration.
#' @param K Number of clusters (a single value or a vector of values)
#' @param replicates The number of iteration for each K.
#' @param model Type of mixture model to use \code{"Poisson"} or 
#' \code{"normal"}. By default, it is the normal.
#' @param GaussianModel Type of \code{GaussianModel} to be used for the 
#' Normal mixture model only. This parameters
#' is set to \code{"Gaussian_pk_Lk_Ck"} by default and doesn't have to 
#' be changed except if an error message proposed
#' to try another model like \code{"Gaussian_pk_Lk_Bk"}.
#' @param transformation The transformation type to be used. By default, 
#' it is the "arcsin" one.
#' @param normFactors The type of estimator to be used to normalize for 
#' differences in library size.
#' By default, it is the "TMM" one.
#' @param merge \code{"union"} or \code{"intersection"}
#' @param clustermq boolean. Does the computation need to be executed 
#' on a distant server. A configuration file and a network connection are 
#' needed to use this option.
#' @param meanFilterCutoff a cutoff to filter a gene with a mean expression lower than this value.
#' (only for RNAseq data, set to NULL for others).
#' @param scale Boolean. If TRUE scale all variables tounit variance. 
#' @param silent if TRUE, coseq run silently (without any console print or 
#' message)
#' @param cmd if TRUE, print steps of the analysis. Used inside the coseq 
#' module in the shiny interface.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: 
#' a tool to process multifactorial RNAseq experiments from quality controls 
#' to co-expression analysis through differential analysis based on contrasts 
#' inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runCoExpression
#' @seealso \code{\link[coseq]{coseq}}
#' @rdname runCoExpression
#' @aliases runCoExpression
#' @importFrom dplyr filter
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"), omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, 
#' formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |> 
#'  runDiffAnalysis(SE.name = "protetest",  
#'                  method = "limmalmFit", contrastList = contrastList) |> 
#'  runCoExpression(SE.name = "protetest", K = 2:5, replicates = 5, 
#'                  contrastNames = contrastList$contrastName, merge = "union")
#'  
#'  getCoexpSettings(MAE[["protetest"]])
#'  
#'  clusters <- getCoseqClusters(MAE[["protetest"]])
#'  
setMethod(f = "runCoExpression",
          signature = "RflomicsSE",
          definition = function(object,
                                K=2:20, 
                                replicates=5, 
                                contrastNames = NULL, 
                                merge="union",
                                model = "Normal",
                                GaussianModel = NULL, 
                                transformation = NULL, 
                                normFactors = NULL, 
                                clustermq = FALSE,
                                meanFilterCutoff = NULL, 
                                scale = NULL,
                                silent = TRUE, 
                                cmd = FALSE){
            
            if (is.null(getDEMatrix(object))) {
              stop("Please run a differential analysis. 
                       runCoExpression uses these results.")
            }
            
            validContrasts <- getValidContrasts(object)
            if(is.null(validContrasts) || nrow(validContrasts) == 0){
              validContrasts <- getSelectedContrasts(object)
              
              if(is.null(validContrasts) || nrow(validContrasts) == 0)
                stop("No defined contrasts")
            }
              
            
            if (is.null(contrastNames)){
              contrastNames <- validContrasts$contrastName
            }
            else{
              contrastNames <- intersect(contrastNames, validContrasts$contrastName)
              if(length(contrastNames) == 0)
                stop("No defined contrasts")
            }
            
            Groups <- getDesignMat(object)
            
            CoExpAnal <- list()
            
            CoExpAnal[["setting"]][["method"]]           <- "coseq"
            CoExpAnal[["setting"]][["gene.list.names"]]  <- contrastNames
            names(CoExpAnal[["setting"]][["gene.list.names"]])  <- contrastNames
            CoExpAnal[["setting"]][["merge.type"]]       <- merge
            CoExpAnal[["setting"]][["replicates.nb"]]    <- replicates
            CoExpAnal[["setting"]][["K.range"]]          <- K
            CoExpAnal[["setting"]][["scale"]]            <- scale
            
            geneList <- getDEList(object = object, 
                                  contrasts = contrastNames,
                                  operation = merge)
            
            # set default parameters based on data type
            param.list <- list("model" = model)
            
            switch(getOmicsTypes(object),
                   
                   "RNAseq" = {
                     counts <- assay(object)[geneList,]
                     
                     param.list[["transformation"]]   <- ifelse(is.null(transformation), "arcsin", transformation)
                     param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "TMM", normFactors)
                     param.list[["meanFilterCutoff"]] <- ifelse(is.null(meanFilterCutoff), 50, meanFilterCutoff)
                     param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Ck", GaussianModel)
                     
                   },
                   "proteomics" = {
                     object2 <- .checkTransNorm(object)
                     counts <- assay(object2)[geneList,]
                     
                     # Print the selected GaussianModel
                     if (cmd) {
                       message("[RFLOMICS] # Use ", GaussianModel)
                       message("[RFLOMICS] # Scale each protein (center = TRUE, scale = TRUE)")
                     } 
                     CoExpAnal[["transformation.prot"]] <- "scaleProt"
                     counts[] <- t(apply(counts,1,function(x){scale(x, center = TRUE, scale = TRUE) }))
                     
                     # param
                     param.list[["transformation"]]   <- ifelse(is.null(transformation), "none", transformation)
                     param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "none", normFactors)
                     param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Bk", GaussianModel)
                   },
                   "metabolomics" = {
                     object2 <- .checkTransNorm(object)
                     counts <-  assay(object2)[geneList,]
                     
                     # Print the selected GaussianModel
                     if (cmd) {
                       message("[RFLOMICS] # Use ", GaussianModel)
                       message("[RFLOMICS] # Scale each metabolite (center = TRUE, scale = TRUE)")
                     } 
                     CoExpAnal[["transformation.metabo"]] <- "scaleMetabo"
                     counts[] <- t(apply(counts,1,function(x){ scale(x, center = TRUE, scale = TRUE) }))
                     
                     # param
                     param.list[["transformation"]] <- ifelse(is.null(transformation), "none", transformation)
                     param.list[["normFactors"]]    <- ifelse(is.null(normFactors), "none", normFactors)
                     param.list[["GaussianModel"]]  <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Bk", GaussianModel)
                   }
            )
            
            CoExpAnal[["setting"]] <- c(CoExpAnal[["setting"]], param.list)
            
            # run coseq : on local machine or remote cluster
            
            if (cmd) message("[RFLOMICS] #     => coseq... ")
            
            
            counts <- counts[, match(rownames(Groups), colnames(counts))]
            if (!identical(colnames(counts), rownames(Groups), attrib.as.set = FALSE)) {
              stop("colnames counts and rownames conds don't match!")
            }
            
            coseq.res.list <- list()
            
            coseq.res.list <- switch(as.character(clustermq),
                                     `FALSE` = {
                                       .tryRflomics(
                                         .runCoseqLocal(counts, 
                                                        conds = Groups$groups,
                                                        K = K, 
                                                        replicates = replicates, 
                                                        param.list = param.list,
                                                        silent = silent,
                                                        cmd = cmd))
                                     },
                                     `TRUE` = {
                                       .tryRflomics(
                                         .runCoseqClustermq(counts, 
                                                            conds = Groups$groups,
                                                            K = K, 
                                                            replicates = replicates, 
                                                            param.list = param.list,
                                                            silent = silent, 
                                                            cmd = cmd))
                                       
                                     })
            
            # If coseq could run (no problem with SSH connexion 
            # in case of clustermq=TRUE)
            
            if(! is.null(coseq.res.list$value)){
              
              CoExpAnal <- c(CoExpAnal, coseq.res.list$value)
            }
            else{
              CoExpAnal[["results"]] <- FALSE
              CoExpAnal[["stats"]]   <- NULL
              CoExpAnal[["error"]]   <- coseq.res.list$error
            }
            
            
            object <- 
              setElementToMetadata(object, 
                                   name = "CoExpAnal", 
                                   content = CoExpAnal)
            return(object)
          })

#' @rdname runCoExpression
#' @aliases runCoExpression
#' @title runCoExpression
#' @exportMethod runCoExpression
setMethod(f          = "runCoExpression",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, K=2:20, 
                                replicates=5, contrastNames, merge="union",
                                model = "Normal", GaussianModel = NULL, 
                                transformation = NULL, normFactors = NULL, 
                                clustermq = FALSE,
                                meanFilterCutoff = NULL, scale = NULL, 
                                silent = TRUE, cmd = FALSE){
            
            
            if (!SE.name %in% names(object)) {
              stop(SE.name, " is not part of ", object)
            }
            
            object[[SE.name]] <- runCoExpression(object = object[[SE.name]],
                                                 K = K,
                                                 replicates = replicates,
                                                 contrastNames = contrastNames, 
                                                 merge = merge, 
                                                 model = model,
                                                 GaussianModel = GaussianModel,
                                                 transformation = transformation,
                                                 normFactors = normFactors,
                                                 clustermq = clustermq,
                                                 meanFilterCutoff = meanFilterCutoff,
                                                 scale = scale,
                                                 silent = silent,
                                                 cmd = cmd)
            
            return(object)
            
          })

#' @title plotCoExpression
#' 
#' @param object An object of class \link{RflomicsSE}
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @importFrom coseq plot
#' @importFrom ggplot2 ggplot aes geom_text geom_boxplot ylim xlab
#' @exportMethod plotCoExpression
#' @rdname plotCoExpression
#' @aliases plotCoExpression
#' @examples
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE(projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"), omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                     typeContrast = "averaged", 
#'                                     returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |> 
#'  runDiffAnalysis(SE.name = "protetest",  method = "limmalmFit", contrastList = contrastList) |> 
#'  runCoExpression(SE.name = "protetest", K = 2:10, replicates = 5, 
#'                  contrastNames = contrastList$contrastName, merge = "union") 
#'  
#'  plots.coexp <- plotCoExpression(MAE[["protetest"]])
#'  plots.coexp[[6]]
#'  
setMethod(f="plotCoExpression",
          signature="RflomicsSE",
          definition <- function(object){
            
            if(is.null(metadata(object)$CoExpAnal) || 
               length( metadata(object)$CoExpAnal) == 0) {
              stop("No co-expression results!")
            }
            CoExpAnal <- metadata(object)$CoExpAnal
            
            Groups <- getDesignMat(object)
            
            coseq.res     <- CoExpAnal[["coseqResults"]]
            ICL.list      <- CoExpAnal[["plots"]][["ICL"]] 
            logLike.list  <- CoExpAnal[["plots"]][["logLike"]]
            K             <- CoExpAnal[["K.range"]]
            
            #### Plots
            ### plot ICL
            limits.vec <- 
              quantile(ICL.list[["ICL.tab"]]$ICL, c(0.1, 0.9), na.rm = TRUE)
            ICL.p <- ggplot(data = ICL.list[["ICL.tab"]]) +
              geom_boxplot(aes(x = as.factor(K), y = ICL, group = K), na.rm = TRUE) +
              scale_y_continuous(limits = limits.vec) +
              geom_text(data = ICL.list[["ICL.n"]], 
                        aes(x = seq_len(length(K)), 
                            y = max(limits.vec),
                            label = paste0("n=", n)), 
                        col = 'red', size = 4) +
              # ylim(min(ICL.list[["ICL.vec"]], na.rm = TRUE),
              #      max(ICL.list[["ICL.vec"]], na.rm = TRUE)) +
              xlab("K")
            
            ### plot logLike
            logLike.p <- ggplot(data = logLike.list[["logLike.tab"]]) +
              geom_boxplot(aes(x = as.factor(K), y = logLike, group = K)) +
              xlab("K") +
              geom_text(data = logLike.list[["logLike.n"]], 
                        aes(x = seq_len(length(K)), 
                            y = max(logLike.list[["logLike.vec"]], na.rm = TRUE),
                            label = paste0("n=", n)), 
                        col = 'red', size = 4)
            
            ### coseq plots
            plot.coseq.res <- plot(coseq.res, 
                                   conds = Groups$groups,
                                   collapse_reps = "average",
                                   graphs = c("profiles", 
                                              "boxplots", 
                                              "probapost_boxplots",
                                              "probapost_barplots",
                                              "probapost_histogram"))
            
            
            return(c(plot.coseq.res, list("ICL" = ICL.p, "logLike" = logLike.p)))
          })

#' @rdname plotCoExpression
#' @aliases plotCoExpression
#' @title plotCoExpression
#' @param SE.name The name of the RflomicsSE to access.
#' @exportMethod plotCoExpression
#' @return A list of plots produced by coseq plot function.
setMethod(f          = "plotCoExpression",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            
            plotCoExpression(object[[SE.name]])
          })



#' @title plotCoExpressionProfile
#'
#' @param object An object of class \link{RflomicsSE}
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @param cluster cluster number
#' @param condition Default is group. 
#' @param features Default is NULL.
#' @return a ggplot graph
#' @exportMethod plotCoExpressionProfile
#' @importFrom dplyr filter mutate rename full_join arrange group_by summarise
#' @importFrom reshape2 melt 
#' @importFrom ggplot2 ggplot aes geom_boxplot aes_string 
#' theme element_text xlab ylab ggtitle geom_point geom_line
#' @rdname plotCoExpressionProfile
#' @aliases plotCoExpressionProfile
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"), omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]],
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |> 
#'  runDiffAnalysis(SE.name = "protetest",  method = "limmalmFit", contrastList = contrastList) |> 
#'  runCoExpression(SE.name = "protetest", K = 2:10, replicates = 5, 
#'                  contrastNames = contrastList$contrastName, merge = "union") 
#'  
#'  plotCoExpressionProfile(MAE[["protetest"]], cluster = 2)
#'  

setMethod(f = "plotCoExpressionProfile",
          signature = "RflomicsSE",
          definition = function(object, cluster = 1, 
                                condition="groups", 
                                features=NULL){
            
            Groups <- getDesignMat(object)
            
            coseq.res  <- metadata(object)$CoExpAnal[["coseqResults"]]
            assays.data <- filter(as.data.frame(coseq.res@assays@data[[1]]), 
                                  get(paste0("Cluster_",cluster)) > 0.8)
            
            y_profiles.gg <- coseq.res@y_profiles[rownames(assays.data),] %>% 
              data.frame() %>% 
              mutate(observations=rownames(.)) %>% 
              melt(id="observations", value.name = "y_profiles") %>%  
              rename(samples = variable) %>%
              full_join(Groups , by = "samples")
            
            y_profiles.gg <- arrange(y_profiles.gg, get(condition))
            y_profiles.gg$groups <- factor(y_profiles.gg$groups, 
                                           levels = unique(y_profiles.gg$groups))
            
            p <- ggplot(data = y_profiles.gg, 
                        aes(x = groups, y = y_profiles)) +
              geom_boxplot(aes_string(fill = condition), 
                           outlier.size = 0.3) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
              xlab("Conditions") + ylab("Expression profiles") +
              ggtitle(paste0("Cluster: ",cluster, "; nb_observations: ",
                             dim(assays.data)[1]))
            
            if(!is.null(features)){
              
              df <- filter(y_profiles.gg, observations == features) %>% 
                group_by(groups) %>% 
                summarise(mean.y_profiles=mean(y_profiles))
              p <- p + 
                geom_point(data = df, 
                           aes(x = groups, y = mean.y_profiles), 
                           color = "red", size = 2) +
                geom_line( data = df, 
                           aes(x = groups, y = mean.y_profiles), 
                           color = "red", group = 1) +
                ggtitle(paste0("Cluster: ",cluster, 
                               "; nb_observations  ", 
                               dim(assays.data)[1], "; red: ", features))
            }                  
            return(p)
          })

#' @rdname plotCoExpressionProfile
#' @title plotCoExpressionProfile
#' @exportMethod plotCoExpressionProfile
#' @aliases plotCoExpressionProfile
setMethod(f          = "plotCoExpressionProfile",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name, cluster = 1, 
                                condition="groups", features=NULL){
            
            plotCoExpressionProfile(object = object[[SE.name]],
                                    cluster = cluster,
                                    condition = condition,
                                    features = observation)
            
          })



#' @title plotCoseqContrasts
#' This function describes the composition of clusters according to the 
#' contrast to which the gene belongs
#' @param object An object of class \link{RflomicsSE}
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @return a ggplot graph
#' @export
#' @exportMethod plotCoseqContrasts
#' @importFrom dplyr filter select starts_with group_by count  ungroup mutate 
#' distinct across left_join
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer 
#' @importFrom ggplot2 ggplot aes  geom_bar  coord_flip labs  
#' scale_fill_discrete theme
#' @rdname plotCoseqContrasts
#' @aliases plotCoseqContrasts
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"),
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"), omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]],
#'                                      typeContrast = "averaged",
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |> 
#'  runDiffAnalysis(SE.name = "protetest",  method = "limmalmFit", contrastList = contrastList) |> 
#'  runCoExpression(SE.name = "protetest", K = 2:10, replicates = 5, 
#'                  contrastNames = contrastList$contrastName, merge = "union") 
#'  
#'  plotCoseqContrasts(MAE[["protetest"]])


setMethod(f = "plotCoseqContrasts",
          signature = "RflomicsSE",
          definition =
            function(object){
              
              # Get Selected contrasts for coexpression
              H <- getCoexpSettings(object)$gene.list.names
              
              if(is.null(H) || length(H) < 2)
                return(NULL)
              
              CoExpAnal <- getAnalysis(object, name = "CoExpAnal")
              
              # Gene's repartition by clusters
              coseq.res  <-
                CoExpAnal[["coseqResults"]]
              genesByclusters <-
                as.data.frame(ifelse(coseq.res@assays@data[[1]] > 0.8, 1, 0))
              genesByclusters <- rownames_to_column(genesByclusters, var = "DEF") 
              
              # Gene's repartition by Contrasts
              genesByContrasts <-
                getDEMatrix(object) |>
                select("DEF",as.vector(H))
              
              # Pivot tab.clusters
              genesByclusters.piv <- pivot_longer(
                data = genesByclusters,
                cols = seq(2, (dim(genesByclusters)[2])),
                names_to = "C",
                values_to = "does.belong") |>
                filter(does.belong == 1) |>
                select(-does.belong)
              
              # Merge the table of Cluster and Contrast
              tab <- left_join(genesByclusters.piv, 
                               as.data.frame(genesByContrasts), by = "DEF") |>
                select(-DEF) |>
                group_by(C) |>
                mutate(sum = rowSums(across(names(genesByContrasts)[-1])))
              
              # Summarize repartition for specific genes
              tab.spe <- filter(tab, sum == 1) |>
                select(-sum) |>
                pivot_longer(names(genesByContrasts)[-1], names_to = "H") |>
                filter(value == 1) |> 
                group_by(C, H) |>
                count()
              
              # Summarize repartition for genes common to at least 
              # 2 contrasts
              tab.com <- filter(tab, sum > 1) |>
                select(-sum) |>
                pivot_longer(names(genesByContrasts)[-1], names_to = "H") |>
                filter(value == 1) |> 
                group_by(C, H) |> 
                count() |> ungroup() |>
                mutate(H = "common") |> distinct()
              
              # Bind table and add prop
              tab <- rbind(tab.com,tab.spe) |> 
                group_by(C) |> 
                mutate(prop=(n/sum(n))*100)
              
              p <-  ggplot(tab,  aes(x = C, y = prop, fill = H)) + 
                geom_bar(stat ="identity") +
                coord_flip() +
                labs(x = "", y = paste0("Proportion of ",
                                        .omicsDic(object)$variableNamegenes)) +
                scale_fill_discrete(name = "",
                                    breaks = c("common", as.vector(H)),
                                    labels = c("commons to at least 2 contrasts", 
                                               paste0("specifics to ", names(H)))) +
                theme(legend.position="bottom",legend.direction = "vertical")
              
              return(p)
            })

#' @rdname plotCoseqContrasts
#' @param SE.name name of the RflomicsSE to access
#' @title plotCoseqContrasts
#' @exportMethod plotCoseqContrasts
#' @aliases plotCoseqContrasts
setMethod(f          = "plotCoseqContrasts",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            
            plotCoseqContrasts(object = object[[SE.name]])
            
          })


# ---- Get Co-exp setting ----
#' @exportMethod getCoexpSettings
#' @rdname Rflomics-accessors
setMethod(f          = "getCoexpSettings",
          signature  = "RflomicsSE",
          
          definition = function(object){
            return(object@metadata$CoExpAnal$setting)   
          })

#' @rdname Rflomics-accessors
#' @exportMethod getCoexpSettings
setMethod(f          = "getCoexpSettings",
          signature  = "RflomicsMAE",
          definition = function(object, SE.name){
            getCoexpSettings(object = object[[SE.name]])
          })

#' @title get members of a cluster
#'
#' @param object a RflomicsSE, produced by rflomics
#' @param clusterName name of the cluster
#' @return The list of entities inside this cluster.
#' @exportMethod getClusterEntities
#' @importFrom coseq clusters
#' @rdname getClusterEntities
#' @aliases getClusterEntities
#' @examples
#' 
#' # Set the data path
#' datPath <- paste0(system.file(package = "RFLOMICS"), "/ExamplesFiles/ecoseed/")
#' 
#' # Load the experimental design
#' ExpDesign <- readExpDesign(file = paste0(datPath, "condition.txt"))
#' 
#' # Set factor name, levels, level of reference and type
#' facRef <- data.frame(factorName   = c("Repeat", "temperature" , "imbibition"), 
#'                      factorRef = c("rep1", "Low", "DS"), 
#'                      factorType = c("batch",  "Bio", "Bio"), 
#'                      factorLevels = c("rep1,rep2,rep3","Low,Medium,Elevated","DS,EI,LI"))
#' 
#' # Load the omics data
#' omicsData <- list( 
#' readOmicsData(file = paste0(datPath, "proteome_ecoseed.txt")))
#'  
#' # Instantiate an object of class RflomicsMAE 
#' MAE <- createRflomicsMAE( projectName = "Tests", omicsData   = omicsData,
#'                          omicsNames  = c("protetest"), omicsTypes  = c("proteomics"),
#'                          ExpDesign   = ExpDesign, factorRef   = facRef)
#' names(MAE) <- c("protetest")
#' 
#' # Set the statistical model and contrasts to test
#' formulae <- generateModelFormulae(MAE)
#' MAE <- setModelFormula(MAE, formulae[[1]])  
#' 
#' # Get the contrasts List and choose the first 3 contrasts of type averaged
#' contrastList <- getPossibleContrasts(MAE, formula = formulae[[1]], 
#'                                      typeContrast = "averaged", 
#'                                      returnTable = TRUE)[c(1, 2, 3),]
#' 
#' # Run the data preprocessing and perform the differential analysis
#'  MAE <- MAE |>  runTransformData(SE.name = "protetest",  transformMethod = "log2") |>
#'  runNormalization(SE.name = "protetest", normMethod = "median")    |> 
#'  runDiffAnalysis(SE.name = "protetest",  method = "limmalmFit", contrastList = contrastList) |> 
#'  runCoExpression(SE.name = "protetest", K = 2:10, replicates = 5, 
#'                  contrastNames = contrastList$contrastName, merge = "union") 
#'  
#'  getClusterEntities(MAE[["protetest"]], clusterName = "cluster_1")
setMethod(
  f          = "getClusterEntities",
  signature  = "RflomicsSE",
  definition = function(object, clusterName) {
    
    clusterName <- gsub("cluster[._]", "", clusterName)
    res <- object@metadata$CoExpAnal$coseqResults
    
    if (!is.null(res)) {
      clList <- clusters(res)
      return(names(which(clList == clusterName)))
    } else {
      return(NULL)
    }
    
  }
)

#' @title get members of a cluster
#'
#' @param object an oblject of class RflomicsMAE, produced by rflomics
#' @param name name of the cluster
#' @param SE.name the name of the data to fetch in the object if the object 
#' is a RflomicsMAE
#' @return The list of entities inside this cluster.
#' @exportMethod getClusterEntities
#' @importFrom coseq clusters
#' @rdname getClusterEntities
#' @aliases getClusterEntities
setMethod(
  f          = "getClusterEntities",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name, clusterName) {
    getClusterEntities(object = object[[SE.name]], 
                       clusterName = clusterName)
  }
)

#' @title get cluster
#'
#' @param object a RflomicsSE, produced by rflomics
#' @return all clusters
#' @exportMethod getCoseqClusters
#' @importFrom coseq clusters
#' @rdname getCoseqClusters
#' @aliases getCoseqClusters
setMethod(
  f          = "getCoseqClusters",
  signature  = "RflomicsSE",
  definition = function(object) {
    
    res <- object@metadata$CoExpAnal$coseqResults
    
    if (!is.null(res)) {
      return(clusters(res))
    } else {
      return(NULL)
    }
  }
)

#' @title get members of a cluster
#'
#' @param object an oblject of class RflomicsMAE, produced by rflomics
#' @param SE.name the name of the data to fetch in the object if the object is 
#' a RflomicsMAE
#' @return all clusters
#' @exportMethod getCoseqClusters
#' @importFrom coseq clusters
#' @rdname getCoseqClusters
#' @aliases getCoseqClusters
setMethod(
  f          = "getCoseqClusters",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name) {
    getCoseqClusters(object = object[[SE.name]])
  }
)


# ---- getCoExpAnalysesSummary ----
#' @title getCoExpAnalysesSummary
#' @description TODO
#' @param object An object of class \link{RflomicsMAE}. 
#' It is expected the SE object is produced by rflomics previous analyses, 
#' as it relies on their results.
#' @param omicNames the name of the experiment to summarize.
#' @param ... Not in use at the moment
#' @return A ggplot object
#' @importFrom dplyr filter mutate rename full_join group_by
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme element_text 
#' facet_grid labs vars
#' @exportMethod getCoExpAnalysesSummary
#' @rdname getCoExpAnalysesSummary
#' @aliases getCoExpAnalysesSummary
setMethod(
  f = "getCoExpAnalysesSummary",
  signature = "RflomicsMAE",
  definition = function(object, omicNames = NULL,
                        ...) {
    mean.y_profiles.list <- list()
    
    if (is.null(omicNames))
      omicNames <- getAnalyzedDatasetNames(object, "CoExpAnal")
    
    for (data in omicNames) {
      
      Groups     <- getDesignMat(object[[data]])
      cluster.nb <-
        object[[data]]@metadata$CoExpAnal$cluster.nb
      coseq.res  <-
        object[[data]]@metadata$CoExpAnal[["coseqResults"]]
      
      mean.y_profiles.list[[data]] <-
        lapply(seq_len(cluster.nb), function(cluster) {
          assays.data <- filter(as.data.frame(coseq.res@assays@data[[1]]),
                                get(paste0("Cluster_", cluster)) > 0.8)
          
          y_profiles.gg <-
            coseq.res@y_profiles[rownames(assays.data), ] %>%
            data.frame() %>%
            mutate(observations = rownames(.)) %>%
            melt(id = "observations", value.name = "y_profiles") %>%
            rename(samples = variable) %>%
            full_join(Groups , by = "samples")
          
          y_profiles.gg %>% group_by(groups) %>%
            summarise(mean = mean(y_profiles)) %>%
            mutate(cluster = paste0("cluster.", cluster))
          
        }) %>% reduce(rbind) %>% mutate(dataset = data)
      
    } %>% reduce(rbind)
    
    mean.y_profiles.gg <- reduce(mean.y_profiles.list, rbind)
    
    if (nrow(mean.y_profiles.gg) == 0)
      return(NULL)
    
    p <- ggplot(data = mean.y_profiles.gg, aes(x = groups, y = mean, group = 1)) +
      geom_line(aes(color = as.factor(cluster))) + 
      geom_point(aes(color = as.factor(cluster))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      facet_grid(cols = vars(cluster), rows = vars(dataset)) +
      labs(x = "Conditions", y = "Expression profiles mean")
    
    return(p)
    
  }
)
