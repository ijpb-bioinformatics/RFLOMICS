#' @title FlomicsExperiment Constructor
#' @description
#' Constructor method for the [\code{\link{FlomicsExperiment-class}}] Class.
#' Create an object of class \code{\link{FlomicsExperiment}} from abundance (rows=features, columns=design)
#' @param abundances numeric matrix or data.frame of features abundance. The names of the columns
#' must give the design experiments using the \code{_} separator (Ex: Fac1_rep1, Fac1_rep2) between
#' factor modalities.
#' @param QCmat An optionnal matrix giving statistics relatives to the quality check (QC) of
#' the experiment. (Ex: For RNAseq: library preparation, sequencing and reads alignments statistics.)
#' @return An object of class [\code{\link{FlomicsExperiment-class}}]
#' @name FlomicsExperiment-Constructor
#' @rdname FlomicsExperiment-Constructor
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
FlomicsExperiment <- function(abundances, QCmat=NULL, ExperimentType="RNAseq", ...)
{
  print("initialize from count matrix and QC matrix if it does exist")

  # Check the class of the count object
  abundances<-as.matrix(abundances)

  # Get the sample names which contain the design
  samples_name <- dimnames(abundances)[[2]]

  # Get the number of design factor and the factors from the names of the matrix count
  design <- GetDesignFromNames(samples_name)

  # Set a list of factor modalities
  dF.List <- lapply(1:design$nb_dFac, function(i){
    as.factor(design$tmpDesign[[i]])
  })
  names(dF.List) <- names(design$tmpDesign)

  # Check if the Quality Check matrix does exist and if true test if the sample names are
  # the same in the two input matrix

  if(! is.null(QCmat)){
  # Check that the names of the factor are the same in the two matrix
    try(if( sum(as.numeric(samples_name %in% row.names(QCmat))) != nrow(QCmat)){
      stop("abundances matrix and QC matrix do not have the same design factors")
    })

  # Get the number of QC factors
  nb_qcFac <- dim(QCmat)[2]
  names(nb_qcFac) <- "n_qcFac"

  # Create a vector with the number of design and QC factors
  nb_dFac <- c(design$nb_dFac, nb_qcFac)

  # Get QC variables
  tmpQC<-QCmat
  names(tmpQC)<-paste("qc",names(QCmat),sep="")

  # Create a colData object with the design factors and QC results vector
  colData <- dplyr::bind_cols(design$tmpDesign,tmpQC)

  }
  else{
    # If the QC matrix does not exist, just return the design factors in the colData object
    colData <- design$tmpDesign
    nb_dFac <- design$nb_dFac
  }

  # instanciate the SummarizedExperiment object
   se <- SummarizedExperiment(assays=
                                S4Vectors::SimpleList(abundances=abundances),
                             colData=colData)
  .FlomicsExperiment(se,
                     ExperimentType=ExperimentType,
                     design=ExperimentalDesign(dF.List=dF.List),
                     colDataStruc=nb_dFac,
                     LogFilter=list(),
                     LogInput=list(),
                     LogTransform=list(),
                     listPCA=list("raw"=FactoMineR::PCA(t(log2(abundances+1)), ncp = 5, graph=F)),
                     Normalization=Normalization(ExperimentType=ExperimentType),
                     DiffAnalysis=new("DiffAnalysis")
                     )
}



#' @title [\code{\link{ExpDesign-class}}] Class constructor
#'
#' @description
#'
#' @param dF.List A list of factor
#'
#' @return An object of class [\code{\link{ExpDesign-class}}]
#' @name ExpDesign-Constructor
#' @rdname ExpDesign-Constructor
#' @export

ExperimentalDesign <- function(dF.List){

  .ExpDesign(List.Factors=dF.List,
           Factors.Type=vector(),
           Model.formula=vector(),
           Contrasts.List=list())
  }


#' @title [\code{\link{Normalization-class}}] Class constructor
#'
#' @description
#'
#' @param ExperimentType The type of experiment (either 'RNAseq' or 'Proteomic' or 'Metabolomic')
#' @return An object of class [\code{\link{Normalization-class}}]
#' @name Normalization-Constructor
#' @rdname Normalization-Constructor
#' @export
#'
Normalization <- function(ExperimentType){

  method=switch(ExperimentType,
                "RNAseq"="TMM",
                "Proteomic"="ProtNormMeth",
                "Metabolomic"="MetaboNormMeth"
                )
  .Normalization(Method=method,
                 Norm.factors=vector())

}



#' @title RunNormalization
#'
#' @param FlomicsExperiment
#'
#' @return FlomicsExperiment
#' @exportMethod RunNormalization
#'
#' @examples
#'
setMethod(f="RunNormalization",
  signature="FlomicsExperiment",
  definition <- function(object, NormMethod=object@Normalization@Method){

      groups <- object@design@List.Factors[object@design@Factors.Type == "Bio"] %>%
                as.data.frame() %>%
                unite(col="groups", sep="_")

      object@Normalization@Norm.factors = switch(NormMethod,
                                                 "TMM"=TMM.Normalization(assay(object), groups$groups)
                                                 )
      return(object)
  }
)



#' @title [\code{\link{DiffAnalysis-class}}] Class constructor
#'
#' @description
#'
#' @return An object of class [\code{\link{DiffAnalysis-class}}]
#' @name DiffAnalysis-Constructor
#' @rdname DiffAnalysis-Constructor
#' @export
DiffAnalysis <- function(){
  .DiffAnalysis(
    method="glmLRT",
    ListOfDEResults=list(),
    ListOfRawMethodResults=list()
  )}


#' @title Multivariate Quality Check
#'
#' @param FlomicsExperiment an object of Class FlomicsExperiment
#' @param axis The number of PCA axis
#'
#' @exportMethod mvQCdesign
#'
#' @rdname mvQCdesign


setMethod(f="mvQCdesign",
          signature="FlomicsExperiment",
          definition <- function(object, data="norm", axis=5){

            #pseudo_abundances <- log2(assay(object)+1)
            #resPCA <- FactoMineR::PCA(t(pseudo_abundances),ncp = axis,graph=F)
            resPCA <- object@listPCA[[data]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <- object@colDataStruc["n_dFac"]

            bigdf <- list()
            for(i in 1:n_dFac){
              Factor <- object@colData[,i]

              df <- list()
              for(j in 1:axis){
                qc = as.vector(resPCA$ind$coord[,j])
                o = order(qc)[order(Factor[order(qc)])]
                col = cc[Factor][o]
                y=qc[o]
                Axis = rep(paste("PCA",j, "\n(Var=",round(resPCA$eig[j,2],1),")",sep=""),length(qc))
                Fac = Factor[o]
                df[[j]] <- data.frame(y,col,Axis,dfac=names(object@colData)[1:n_dFac][i],
                                      Levels=Fac,x=1:length(y))
              }
              bigdf[[i]] <- bind_rows(df)
            }
            big <- dplyr::bind_rows(bigdf)
            out <- by(data = big, INDICES = big$dfac, FUN = function(m) {
              m <- droplevels(m)
              m <- ggplot(m,aes(y=y,x=x,color=Levels))+
                geom_bar(stat = "identity",position = position_dodge(),aes(fill=Levels))+
                facet_grid(as.factor(dfac)~Axis) +
                labs(x = "Samples", y="Coordinates on \n the PCA axis")+
                theme(axis.title.y = element_text(size = 5),
                      axis.title.x=element_text(size = 5),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
            })
            do.call(grid.arrange, out)
          })


#' @title multivariate QC data
#'
#' @param FlomicsExperiment an object of Class FlomicsExperiment
#' @param axis The number of PCA axis
#'
#' @exportMethod mvQCdata
#' @rdname mvQCdata

setMethod(f="mvQCdata",
          signature="FlomicsExperiment",
          definition <- function(object,axis=3){

            pseudo_abundances <- log2(assay(object)+1)
            resPCA <- FactoMineR::PCA(t(pseudo_abundances),ncp = axis,graph=F)
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <- object@colDataStruc["n_dFac"]
            n_qcFac <- object@colDataStruc["n_qcFac"]
            var <- names(object@colData[,(n_dFac+1):(n_dFac+n_qcFac)])

            corA=list()
            VarAxis = list()
            for(i in 1:axis){
            corA[[i]] <- cor(resPCA$ind$coord[,i],
                             as.data.frame(object@colData[,(n_dFac+1):(n_dFac+n_qcFac)]),
                             method="spearman")
            VarAxis[[i]] <- paste("\n(Var=",round(resPCA$eig[i,2],1),")",sep="")
            }

            df = data.frame("QCparam"=rep(var,axis),
                            "Spearman"= unlist(corA),
                            "Axis"=rep(paste(rep("Axis",axis),1:axis,VarAxis,sep=""), each=n_qcFac))

            ggplot(df,aes(x=Axis, y=abs(Spearman),fill=QCparam))+
              geom_bar(stat="identity",position=position_dodge(),width=0.7)+ylim(0,1)+
              labs(x = "Axis number", y="Cor(Coord_dFactor_PCA,QCparam)")
          })



#' @title boxplotQCnorm
#'
#' @param FlomicsExperiment an object of Class FlomicsExperiment
#'
#' @exportMethod boxplotQCnorm
#' @rdname boxplotQCnorm

setMethod(f= "boxplotQCnorm",
          signature = "FlomicsExperiment",
          definition <- function(object){

            # this function generate boxplot (abandance distribution) from raw data and normalized data

            #col <- colorPlot(object@design, object@colData, condition="samples")
            sample_names <- object@design@List.Factors %>% as.data.frame() %>%
                            unite(., col="samples", sep="_")
            groups  <- object@design@List.Factors[object@design@Factors.Type == "Bio"] %>% as.data.frame() %>%
                       unite(col="groups", sep="_", remove = FALSE) %>% mutate(samples=sample_names$samples)

            # raw data
            pseudo_raw <- log2(assay(object)+1) %>%  data.table::melt() %>% mutate(TAG = "1.Raw data")

            # normalized data
            pseudo_norm  <- log2(scale(assay(object),center=FALSE,scale=object@Normalization@Norm.factors$norm.factors)+1) %>%
                                    data.table::melt() %>% mutate(TAG = "2.Normalized data")

            pseudo_tmp <- rbind(pseudo_raw, pseudo_norm)
            # merge  raw data & normalized data

            colnames(pseudo_tmp) <- c("features", "samples", "counts", "TAG")

            pseudo <- full_join(pseudo_tmp, groups, by="samples")

            pseudo$samples <- factor(pseudo$samples, levels = unique(pseudo$samples))

            # boxplot
            ggplot(pseudo, aes(x=samples, y=counts)) +
              geom_boxplot(aes(fill=groups)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              facet_grid(~TAG)
              #scale_fill_manual(values=col)

          }
)


#' @title plotPCAnorm
#'
#' @param FlomicsExperiment
#' @param condition
#' @param color color palette
#'
#' @return
#' @exportMethod plotPCAnorm
#'
#' @examples
setMethod(f= "plotPCAnorm",
          signature = "FlomicsExperiment",
          definition <- function(object, data="norm", PCs=c(1,2), condition="groups"){

            #
            PC1 <- paste("Dim.",PCs[1], sep="")
            PC2 <- paste("Dim.",PCs[2], sep="")

            #col <- colorPlot(object@design, object@colData, condition=condition)

            #sample_names <- row.names(object@listPCA[[data]]$ind$coord)

            groups  <- object@design@List.Factors[object@design@Factors.Type == "Bio"] %>% as.data.frame() %>%
                       unite(col="groups", sep="_")
            factors <- object@design@List.Factors %>% as.data.frame() %>%
                       unite(., col="samples", sep="_", remove = FALSE) %>% mutate(groups = groups$groups)

            score_tmp <- object@listPCA[[data]]$ind$coord[, PCs] %>% as.data.frame() %>% mutate(samples=row.names(.))

            score  <- full_join(score_tmp, factors, by="samples")

            ggplot(score, aes_string(x=PC1, y=PC2, color=condition))  +
              geom_point(size=3) +
              geom_text(aes(label=samples), size=3, vjust = 0) +
              xlab(paste(PC1, " (",round(object@listPCA[[data]]$eig[PCs,2][1], digit=3),"%)", sep="")) +
              ylab(paste(PC2, " (",round(object@listPCA[[data]]$eig[PCs,2][2], digit=3),"%)", sep="")) +
              geom_hline(yintercept=0, linetype="dashed", color = "red") +
              geom_vline(xintercept=0, linetype="dashed", color = "red") +
              theme(strip.text.x = element_text(size=8, face="bold.italic"),
                    strip.text.y = element_text(size=8, face="bold.italic")) #+
              #scale_color_manual(values=col$colors)
            })



#' @title barplotPCAnorm
#'
#' @param FlomicsExperiment
#' @param condition
#' @param colors color palette
#'
#' @return
#' @exportMethod barplotPCAnorm
#' 
#'
#' @examples
setMethod(f= "barplotPCAnorm",
          signature = "FlomicsExperiment",
          definition <- function(object, condition="samples"){

            col <- colorPlot(object@design, object@colData, condition=condition)

            score_raw  <- object@listPCA$raw$ind$coord  %>% data.table::melt %>%
                          mutate(tag="1.Unnormalised data")
            score_norm <- object@listPCA$norm$ind$coord %>% data.table::melt %>%
                          mutate(tag=paste("2.Normalised data : ", object@Normalization@Method,  sep=""))

            score <- rbind(score_raw, score_norm)
            colnames(score) <- c("samples", "PCs", "value", "tag")

            ggplot(data=score, aes(x=PCs, y=value, fill=samples)) +
              geom_bar(stat="identity", position=position_dodge(), color="black") +
              facet_grid(tag~PCs, scale ="free", space = "free") + scale_fill_manual(values=col)
          })



#' FilterLowAbundance
#'
#' @param FlomicsExperiment
#' @param threshold
#'
#' @return FlomicsExperiment
#' @exportMethod FilterLowAbundance
#'
#' @examples
setMethod(f= "FilterLowAbundance",
          signature = "FlomicsExperiment",
          definition <- function(object, threshold){


            feature_0 <- object[rowSums(assay(object)) <= threshold, ]@NAMES
            
            object    <- object[rowSums(assay(object)) > threshold, ]
            
            BioFact   <- names(object@design@List.Factors[object@design@Factors.Type == "Bio"])
            
            #Replicat  <- levels(object@design@List.Factors[object@design@Factors.Type != "Bio"][[1]])

            object@LogFilter[["feature_0"]] <- feature_0
            #object@LogFilter[["current"]]   <- data.frame(number   =c(dim(assay(object)), length(BioFact), length(Replicat)), 
            #                                              row.names=c("Features", "Samples", "Bio Factors", "Replicats"))
            object@LogFilter[["current"]]   <- data.frame(number   =c(dim(assay(object)), length(BioFact)), 
                                                          row.names=c("Features", "Samples", "Bio Factors"))
            
            return(object)
          })
