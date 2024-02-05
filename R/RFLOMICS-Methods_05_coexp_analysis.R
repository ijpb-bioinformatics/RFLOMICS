
################################### CO-EXPRESSION #############################



#' @title runCoExpression
#' @description This is an interface method which performs co-expression/co-abundance analysis
#' of omic-data.
#' @details For now, only the coseq function of the coseq package is used.
#' For RNAseq data, parameters used are those recommended in DiCoExpress workflow (see the reference).
#' This parameters are: \code{model="normal"}, \code{transformation="arcsin"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="TMM"}, \code{meanFilterCutoff = 50}
#' For proteomic or metabolomic, data are scaled by protein or metabolite to group them by expression
#' profiles rather than by expression intensity.
#' After data scaling, recommended parameters (from \code{coseq} developers) for co-expression analysis are:
#' \code{model="normal"}, \code{transformation="none"}, \code{GaussianModel="Gaussian_pk_Lk_Ck"},
#' \code{normFactors="none"},  \code{meanFilterCutoff = NULL}.
#'
#' @return
#' An S4 object of class \link{RflomicsSE}
#' All the results are stored as a named list \code{CoExpAnal} in the metadata slot of a
#' given \code{RflomicsSE} object. Objects are:
#' The runCoExpression method return several results, for \link{coseq} method, objects are:
#' \itemize{
#' \item{\code{model:} }{see model params description}
#' \item{\code{transformation:} }{see transformation params description}
#' \item{\code{normFactors:} }{see normFactors params description}
#' \item{\code{meanFilterCutoff:} }{set to 50 for RNA and to NULL for others}
#' \item{\code{gene.list.names:} }{see nameList in Arguments description}
#' \item{\code{merge.type:} }{see merge params description}
#' \item{\code{coseqResults:} }{the raw results of \code{coseq}}
#' \item{\code{clusters:} }{a List of clusters}
#' \item{\code{cluster.nb:} }{The number of cluster}
#' \item{\code{plots:} }{The plots of \code{coseq} results}
#' }
#' @param object An object of class \link{RflomicsSE}
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @param nameList names of the contrasts from which the DE entities are taken. Can be NULL, in that case every contrasts from the differential analysis is taken into consideration.
#' @param K Number of clusters (a single value or a vector of values)
#' @param replicates The number of iteration for each K.
#' @param model Type of mixture model to use \code{"Poisson"} or \code{"normal"}. By default, it is the normal.
#' @param GaussianModel Type of \code{GaussianModel} to be used for the Normal mixture model only. This parameters
#' is set to \code{"Gaussian_pk_Lk_Ck"} by default and doesn't have to be changed except if an error message proposed
#' to try another model like \code{"Gaussian_pk_Lk_Bk"}.
#' @param transformation The transformation type to be used. By default, it is the "arcsin" one.
#' @param normFactors The type of estimator to be used to normalize for differences in library size.
#' By default, it is the "TMM" one.
#' @param merge \code{"union"} or \code{"intersection"}
#' @param clustermq_arg boolean. Does the computation need to be executed on a distant server?
#' @param silent if TRUE, coseq run silently (without any console print or message)
#' @param cmd if TRUE, print steps of the analysis. Used inside the coseq module in the shiny interface.
#' @references
#' Lambert, I., Paysant-Le Roux, C., Colella, S. et al. DiCoExpress: a tool to process multifactorial RNAseq experiments from quality controls to co-expression analysis through differential analysis based on contrasts inside GLM models. Plant Methods 16, 68 (2020).
#' @exportMethod runCoExpression
#' @seealso \code{\link{coseq::coseq}}
#' @rdname runCoExpression
#' 
methods::setMethod(f = "runCoExpression",
                   signature = "RflomicsSE",
                   definition = function(object,
                                         K=2:20, 
                                         replicates=5, 
                                         nameList = NULL, 
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
                     
                     if (is.null(object@metadata$DiffExpAnal[["mergeDEF"]]))
                       stop("Please run a differential analysis. runCoExpression uses these results.")
                     
                     if (is.null(nameList) && !is.null(getValidContrasts(object)[["tag"]])) 
                       nameList <- getValidContrasts(object)[["tag"]]
                     else if (is.null(nameList) && is.null(getValidContrasts(object)[["tag"]])) 
                       nameList <- colnames(object@metadata$DiffExpAnal[["mergeDEF"]])[-1]
                     
                     Groups <- getDesignMat(object)
                     
                     CoExpAnal <- list()
                     
                     CoExpAnal[["setting"]][["method"]]           <- "coseq"
                     CoExpAnal[["setting"]][["gene.list.names"]]  <- nameList
                     names(CoExpAnal[["setting"]][["gene.list.names"]])  <- dplyr::filter(object@metadata$DiffExpAnal$contrasts, tag %in% nameList)$contrastName
                     CoExpAnal[["setting"]][["merge.type"]]       <- merge
                     CoExpAnal[["setting"]][["replicates.nb"]]    <- replicates
                     CoExpAnal[["setting"]][["K.range"]]          <- K
                     CoExpAnal[["setting"]][["scale"]]            <- scale
                     
                     geneList <- getDEList(object = object, contrasts = nameList, operation = merge)
                     
                     # set default parameters based on data type
                     param.list <- list("model" = model)
                     
                     switch(object@metadata$omicType,
                            
                            "RNAseq" = {
                              counts <- SummarizedExperiment::assay(object)[geneList,]
                              
                              param.list[["transformation"]]   <- ifelse(is.null(transformation), "arcsin", transformation)
                              param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "TMM", normFactors)
                              param.list[["meanFilterCutoff"]] <- ifelse(is.null(meanFilterCutoff), 50, meanFilterCutoff)
                              param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Ck", GaussianModel)
                              
                            },
                            "proteomics" = {
                              object <- checkTransNorm(object)
                              counts <- SummarizedExperiment::assay(object)[geneList,]
                              
                              # Print the selected GaussianModel
                              if (cmd) print(paste("Use ", GaussianModel, sep = ""))
                              if (cmd) print("Scale each protein (center = TRUE, scale = TRUE)")
                              CoExpAnal[["transformation.prot"]] <- "scaleProt"
                              counts[] <- t(apply(counts,1,function(x){scale(x, center = TRUE, scale = TRUE) }))
                              
                              # param
                              param.list[["transformation"]]   <- ifelse(is.null(transformation), "none", transformation)
                              param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "none", normFactors)
                              param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Bk", GaussianModel)
                            },
                            "metabolomics" = {
                              object <- checkTransNorm(object)
                              counts <- SummarizedExperiment::assay(object)[geneList,]
                              
                              # Print the selected GaussianModel
                              if (cmd) print(paste("Use ", GaussianModel, sep= ""))
                              if (cmd) print("Scale each metabolite (center = TRUE,scale = TRUE)")
                              CoExpAnal[["transformation.metabo"]] <- "scaleMetabo"
                              counts[] <- t(apply(counts,1,function(x){ scale(x, center = TRUE, scale = TRUE) }))
                              
                              # param
                              param.list[["transformation"]]   <- ifelse(is.null(transformation), "none", transformation)
                              param.list[["normFactors"]]      <- ifelse(is.null(normFactors), "none", normFactors)
                              param.list[["GaussianModel"]]    <- ifelse(is.null(GaussianModel), "Gaussian_pk_Lk_Bk", GaussianModel)
                            }
                     )
                     
                     CoExpAnal[["setting"]] <- c(CoExpAnal[["setting"]], param.list)
                     
                     # run coseq : on local machine or remote cluster
                     
                     if (cmd) print("#     => coseq... ")
                     
                     
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
                     
                     # If coseq could run (no problem with SSH connexion in case of clustermq=TRUE)
                     
                     if(! is.null(coseq.res.list$value)){
                       
                       CoExpAnal <- c(CoExpAnal, coseq.res.list$value)
                     }
                     else{
                       CoExpAnal[["results"]] <- FALSE
                       CoExpAnal[["stats"]]   <- NULL
                       CoExpAnal[["error"]]   <- coseq.res.list$error
                     }
                     
                     object@metadata$CoExpAnal <- CoExpAnal
                     return(object)
                   })

#' @rdname runCoExpression
#' @title runCoExpression
#' @exportMethod runCoExpression
methods::setMethod(f          = "runCoExpression",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, K=2:20, replicates=5, nameList, merge="union",
                                         model = "Normal", GaussianModel = NULL, 
                                         transformation = NULL, normFactors = NULL, clustermq = FALSE,
                                         meanFilterCutoff = NULL, scale = NULL, silent = TRUE, cmd = FALSE){
                     
                     
                     if (!SE.name %in% names(object)) 
                       stop(paste0(SE.name, " is not part of ", object))
                     
                     object[[SE.name]] <-  runCoExpression(object = object[[SE.name]],
                                                           K = K,
                                                           replicates = replicates,
                                                           nameList = nameList, 
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

# Pour utiliser la fonction repeatable(), "seed"  pourrait être ajouté en paramètre.




#' @title plotCoExpression
#' 
#' @param object An object of class \link{RflomicsSE}
#' @return list plot of ICL, logLike and coseq object with min ICL
#' @importFrom coseq plot
#' @importFrom ggplot2 ggplot geom_boxplot geom_text
#' @exportMethod plotCoExpression
#' @noRd
#' 
methods::setMethod(f="plotCoExpression",
                                        signature="RflomicsSE",
                                        definition <- function(object){
                                          
                                          if(is.null(object@metadata$CoExpAnal) || length(object@metadata$CoExpAnal) == 0) stop("No co-expression results!")
                                          CoExpAnal <- object@metadata$CoExpAnal
                                          
                                          Groups <- getDesignMat(object)
                                          
                                          coseq.res     <- CoExpAnal[["coseqResults"]]
                                          ICL.list      <- CoExpAnal[["plots"]][["ICL"]] 
                                          logLike.list  <- CoExpAnal[["plots"]][["logLike"]]
                                          K             <- CoExpAnal[["K.range"]]
                                          
                                          #### Plots
                                          ### plot ICL
                                          ICL.p   <- ggplot2::ggplot(data = ICL.list[["ICL.tab"]]) +
                                            ggplot2::geom_boxplot(ggplot2::aes(x = as.factor(K), y = ICL, group = K)) +
                                            ggplot2::geom_text(data = ICL.list[["ICL.n"]], ggplot2::aes(x = 1:length(K), y = max(ICL.list[["ICL.vec"]], na.rm = TRUE),
                                                                                                        label = paste0("n=", n)), col = 'red', size = 4) +
                                            ggplot2::ylim(min(ICL.list[["ICL.vec"]], na.rm = TRUE), max(ICL.list[["ICL.vec"]], na.rm = TRUE)) +
                                            ggplot2::xlab("K")
                                          
                                          ### plot logLike
                                          logLike.p   <- ggplot2::ggplot(data = logLike.list[["logLike.tab"]]) +
                                            ggplot2::geom_boxplot(ggplot2::aes(x = as.factor(K), y = logLike, group = K)) +
                                            ggplot2::xlab("K") +
                                            ggplot2::geom_text(data = logLike.list[["logLike.n"]], ggplot2::aes(x = 1:length(K), y = max(logLike.list[["logLike.vec"]], na.rm = TRUE),
                                                                                                                label = paste0("n=", n)), col = 'red', size = 4)
                                          
                                          ### coseq plots
                                          plot.coseq.res <- coseq::plot(coseq.res, conds = Groups$groups, collapse_reps = "average",
                                                                        graphs = c("profiles", "boxplots", "probapost_boxplots",
                                                                                   "probapost_barplots", "probapost_histogram"))
                                          
                                          # CoExpAnal[["plots"]] <- plot.coseq.res
                                          # CoExpAnal[["plots"]][["ICL"]]     <- ICL.p
                                          # CoExpAnal[["plots"]][["logLike"]] <- logLike.p
                                          
                                          return(c(plot.coseq.res, list("ICL" = ICL.p, "logLike" = logLike.p)))
                                        })

#' @rdname plotCoExpression
#' @title plotCoExpression
#' @exportMethod plotCoExpression
methods::setMethod(f          = "plotCoExpression",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, numCluster = 1, condition="groups", observation=NULL){
                     
                     plotCoExpression(object = object[[SE.name]],
                                             numCluster = numCluster,
                                             condition = condition,
                                             observation = observation)
                   })



#' @title plotCoExpressionProfile
#'
#' @param object An object of class \link{RflomicsSE}
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @param numCluster cluster number
#' @param condition 
#' @param observation 
#' @export
#' @exportMethod plotCoExpressionProfile
#' @importFrom dplyr filter mutate rename full_join arrange group_by summarise
#' @importFrom reshape2 melt 
#' @noRd

methods::setMethod(f="plotCoExpressionProfile",
                   signature="RflomicsSE",
                   definition <- function(object, numCluster = 1, condition="groups", observation=NULL){
                     
                     Groups <- getDesignMat(object)
                     
                     coseq.res  <- object@metadata$CoExpAnal[["coseqResults"]]
                     assays.data <- dplyr::filter(as.data.frame(coseq.res@assays@data[[1]]), get(paste0("Cluster_",numCluster)) > 0.8)
                     
                     y_profiles.gg <- coseq.res@y_profiles[rownames(assays.data),] %>% 
                       data.frame() %>% 
                       dplyr::mutate(observations=rownames(.)) %>% 
                       reshape2::melt(id="observations", value.name = "y_profiles") %>%  
                       dplyr::rename(samples = variable) %>%
                       dplyr::full_join(Groups , by = "samples")
                     
                     y_profiles.gg <- dplyr::arrange(y_profiles.gg, get(condition))
                     y_profiles.gg$groups <- factor(y_profiles.gg$groups, levels = unique(y_profiles.gg$groups))
                     
                     p <- ggplot2::ggplot(data = y_profiles.gg, ggplot2::aes(x = groups, y = y_profiles)) +
                       ggplot2::geom_boxplot(ggplot2::aes_string(fill = condition), outlier.size = 0.3) +
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +  
                       ggplot2::xlab("Conditions") + ggplot2::ylab("Expression profiles") +
                       ggplot2::ggtitle(paste0("Cluster: ",numCluster, "; nb_observations : ", dim(assays.data)[1]))
                     
                     if(!is.null(observation)){
                       
                       df <- dplyr::filter(y_profiles.gg, observations == observation) %>% 
                         dplyr::group_by(groups) %>% 
                         dplyr::summarise(mean.y_profiles=mean(y_profiles))
                       p <- p + 
                         ggplot2::geom_point(data = df, ggplot2::aes(x = groups, y = mean.y_profiles), color = "red", size = 2) +
                         ggplot2::geom_line( data = df, ggplot2::aes(x = groups, y = mean.y_profiles), color = "red", group = 1) +
                         ggplot2::ggtitle(paste0("Cluster: ",numCluster, "; nb_observations : ", dim(assays.data)[1], "; red : ", observation))
                     }
                     -                     
                       return(p)
                   })

#' @rdname plotCoExpressionProfile
#' @title plotCoExpressionProfile
#' @exportMethod plotCoExpressionProfile
methods::setMethod(f          = "plotCoExpressionProfile",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name, numCluster = 1, condition="groups", observation=NULL){
                     
                     plotCoExpressionProfile(object = object[[SE.name]],
                                        numCluster = numCluster,
                                        condition = condition,
                                        observation = observation)
                     
                   })



#' @title plotCoseqContrasts
#' This function describes the composition of clusters according to the contrast to which the gene belongs
#' @param object An object of class \link{RflomicsSE}
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @export
#' @exportMethod plotCoseqContrasts
#' @importFrom dplyr filter mutate rename full_join arrange group_by summarise
#' @noRd


methods::setMethod(f = "plotCoseqContrasts",
                   signature = "RflomicsSE",
                   definition <-
                     function(object){
                       
                       # Get Selected contrasts for coexpression
                       H <- getCoexpSetting(object)$gene.list.names
                       
                       # Gene's repartition by clusters
                       coseq.res  <-
                         object@metadata$CoExpAnal[["coseqResults"]]
                       genesByclusters <-
                         as.data.frame(ifelse(coseq.res@assays@data[[1]] > 0.8, 1, 0))
                       genesByclusters <-  tibble::rownames_to_column(genesByclusters,var="DEF") 
                       
                       # Gene's repartition by Contrasts
                       genesByContrasts <-
                         object@metadata$DiffExpAnal$mergeDEF |>
                         dplyr::select("DEF",as.vector(H))
                       
                       # Pivot tab.clusters
                       genesByclusters.piv <- tidyr::pivot_longer(
                         data = genesByclusters,
                         cols = 2:(dim(genesByclusters)[2]),
                         names_to = "C",
                         values_to = "does.belong") |>
                         dplyr::filter(does.belong == 1) |>
                         dplyr::select(-does.belong)
                       
                       # Merge the table of Cluster and Contrast
                       tab <- dplyr::left_join(genesByclusters.piv, as.data.frame(genesByContrasts),by="DEF") |>
                         dplyr::select(-DEF) |>
                         dplyr::group_by(C) |>
                         dplyr::mutate(sum = rowSums(dplyr::across(dplyr::starts_with("H"))))
                       
                       # Summarize repartition for specific genes
                       tab.spe <- dplyr::filter(tab, sum == 1) |>
                         dplyr::select(-sum) |>
                         tidyr::pivot_longer(dplyr::starts_with("H"), names_to = "H") |>
                         dplyr::filter(value == 1) |> 
                         dplyr::group_by(C, H) |>
                         dplyr::count()
                       
                       # Summarize repartition for genes common to at least 2 contrasts
                       tab.com <- dplyr::filter(tab, sum > 1) |>
                         dplyr::select(-sum) |>
                         tidyr::pivot_longer(dplyr::starts_with("H"), names_to = "H") |>
                         dplyr::filter(value == 1) |> 
                         dplyr::group_by(C, H) |> 
                         dplyr::count() |> dplyr::ungroup() |>
                         dplyr::mutate(H= "common") |> dplyr::distinct()
                       
                       # Bind table and add prop
                       tab <- rbind(tab.com,tab.spe) |> group_by(C) |> mutate(prop=(n/sum(n))*100)
                       
                       p <- ggplot2::ggplot(tab, ggplot2::aes(x = C, y = prop, fill = H)) + 
                         ggplot2::geom_bar(stat ="identity") +
                         ggplot2::coord_flip() +
                         ggplot2::labs(x="",y=paste0("Proportion of", RFLOMICS:::omicsDic(object)$variableNamegenes)) +
                         ggplot2::scale_fill_discrete(name="",
                                                      breaks=c("common", as.vector(H)),
                                                      labels=c("commons to at least 2 contrasts", 
                                                               paste0("specifics to ",names(H)))) +
                         ggplot2::theme(legend.position="bottom",legend.direction = "vertical")
                       
                       return(p)
                     })

#' @rdname plotCoseqContrasts
#' @title plotCoseqContrasts
#' @exportMethod plotCoseqContrasts
methods::setMethod(f          = "plotCoseqContrasts",
                   signature  = "RflomicsMAE",
                   definition = function(object){
                     
                     plotCoseqContrasts(object = object[[SE.name]])
                     
                   })


# ---- Get diff setting ----

#' @title Get differential analysis setting parameters
#'
#' @param object of class RflomicsSE
#' @return List of differential analysis setting parametres.
#' @exportMethod getCoexpSetting
#' @rdname getCoexpSetting
#'

methods::setMethod(f          = "getCoexpSetting",
                   signature  = "RflomicsSE",
                   
                   definition = function(object){
                     return(object@metadata$CoExpAnal$setting)   
                   })

#' @rdname getCoexpSetting
#' @title getCoexpSetting
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @exportMethod getCoexpSetting

methods::setMethod(f          = "getCoexpSetting",
                   signature  = "RflomicsMAE",
                   definition = function(object, SE.name){
                     getCoexpSetting(object = object[[SE.name]])
                   })



#' # ---- INTERNAL - get members of a cluster or coseq clusters ----
#' #
#' #' @title get members of a cluster
#' #'
#' #' @param object a RflomicsSE, produced by rflomics
#' #' @param name name of the cluster
#' #' @return The list of entities inside this cluster.
#' #' @noRd
#' #' @importFrom coseq clusters
#' #' @keywords internal
#' 
#' .getCluster <- function(object, clusterName) {
#'   
#'   clusterName <- gsub("cluster[.]", "", clusterName)
#'   res <- object@metadata$CoExpAnal$coseqResults
#'   
#'   if (!is.null(res)) {
#'     clList <- clusters(res)
#'     return(names(clList == clusterName))
#'   } else {
#'     return(NULL)
#'   }
#'   
#' }

#' @title get members of a cluster
#'
#' @param object a RflomicsSE, produced by rflomics
#' @param clusterName name of the cluster
#' @return The list of entities inside this cluster.
#' @exportMethod getClusterEntities
#' @importFrom coseq clusters


methods::setMethod(
  f          = "getClusterEntities",
  signature  = "RflomicsSE",
  definition = function(object, clusterName) {
    
    clusterName <- gsub("cluster[.]", "", clusterName)
    res <- object@metadata$CoExpAnal$coseqResults
    
    if (!is.null(res)) {
      clList <- clusters(res)
      return(names(clList == clusterName))
    } else {
      return(NULL)
    }
    
  }
)

#' @title get members of a cluster
#'
#' @param object an oblject of class RflomicsMAE, produced by rflomics
#' @param name name of the cluster
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @return The list of entities inside this cluster.
#' @exportMethod getClusterEntities
#' @importFrom coseq clusters


methods::setMethod(
  f          = "getClusterEntities",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name, clusterName) {
    getClusterEntities(object = object[[SE.name]], clusterName = clusterName)
  }
)

#' #' @param object a RflomicsSE, produced by rflomics
#' #' @return all clusters
#' #' @noRd
#' #' @importFrom coseq clusters
#' #' @keywords internal
#' 
#' .getCoseqClusters <- function(object) {
#'   
#'   res <- object@metadata$CoExpAnal$coseqResults
#'   
#'   if (!is.null(res)) {
#'     return(clusters(res))
#'   } else {
#'     return(NULL)
#'   }
#'   
#' }


#' @title get cluster
#'
#' @param object a RflomicsSE, produced by rflomics
#' @return all clusters
#' @exportMethod getCoseqClusters
#' @importFrom coseq clusters


methods::setMethod(
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
#' @param SE.name the name of the data to fetch in the object if the object is a RflomicsMAE
#' @return all clusters
#' @exportMethod getCoseqClusters
#' @importFrom coseq clusters


methods::setMethod(
  f          = "getCoseqClusters",
  signature  = "RflomicsMAE",
  definition = function(object, SE.name) {
    getCoseqClusters(object = object[[SE.name]])
  }
)


