#' @title FlomicsExperiment Constructor
#' @description
#' Constructor method for the [\code{\link{FlomicsExperiment-class}}] Class.
#' Create an object of class \code{\link{FlomicsExperiment}} from counts (rows=features, columns=design)
#' @param counts numeric matrix or data.frame of read counts. The names of the columns
#' must give the design experiments using the \code{_} separator (Ex: Fac1_rep1, Fac1_rep2) between
#' factor modalities.
#' @param QCmat An optionnal matrix giving statistics relatives to the quality check (QC) of
#' the RNAseq experiment: library preparation, sequencing and reads alignments statistics.
#' @return An object of class [\code{\link{FlomicsExperiment-class}}]
#' @name FlomicsExperiment-Constructor
#' @rdname FlomicsExperiment-Constructor
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
FlomicsExperiment <- function(counts, QCmat=NULL, ...)
{
  print("initialize from count matrix and QC matrix if it does exist")

  # Check the class of the count object
  counts<-as.matrix(counts)

  # Get the sample names which contain the design
  samples_name <- dimnames(counts)[[2]]

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
      stop("Counts matrix and QC matrix do not have the same design factors")
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
                                S4Vectors::SimpleList(counts=counts),
                             colData=colData)
  .FlomicsExperiment(se, colDataStruc=nb_dFac, design=ExperimentalDesign(dF.List=dF.List))
}



ExperimentalDesign <- function(dF.List){

  .ExpDesign(List.Factors=dF.List,
           Factors.Type=vector(),
           Model.formula=formula(),
           Contrasts.List=list())
  }




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
          definition <- function(object,axis=3){

            pseudo_counts <- log2(assay(object)+1)
            resPCA <- FactoMineR::PCA(t(pseudo_counts),ncp = axis,graph=F)
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
              df[[j]] <- data.frame(y,col,Axis,dfac=names(object@colData)[1:n_dFac][i],Levels=Fac,x=1:length(y))
            }
            bigdf[[i]] <- bind_rows(df)
            }
            big <- dplyr::bind_rows(bigdf)
            ggplot(big,aes(y=y,x=x,color=Levels))+
              geom_bar(stat = "identity",position = position_dodge())+
              facet_grid(as.factor(dfac)~Axis) +
              labs(x = "Factors", y="Coordinates on the PCA axis")
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

            pseudo_counts <- log2(assay(object)+1)
            resPCA <- FactoMineR::PCA(t(pseudo_counts),ncp = axis,graph=F)
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

