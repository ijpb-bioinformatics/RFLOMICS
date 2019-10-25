#' @title [\code{\link{ExpDesign-class}}] Class constructor
#' @description
#' @param dF.List A list of factor
#' @return An object of class [\code{\link{ExpDesign-class}}]
#' @name ExpDesign-Constructor
#' @rdname ExpDesign-Constructor
#' @export
ExperimentalDesign <- function(ExpDesign){

  dF.List <- lapply(1:dim(ExpDesign)[2], function(i){
    as.factor(ExpDesign[[i]])
  })
  names(dF.List) <- names(ExpDesign)
  
  .ExpDesign(List.Factors=dF.List,
           Factors.Type=vector(),
           Model.formula=vector(),
           Contrasts.List=list())
  }



#' @title RunNormalization
#' @param MultiAssayExperiment
#' @param data omic data type
#' @param NormMethod normalisation methode
#' @return MultiAssayExperiment
#' @exportMethod RunNormalization
#' @examples
#' 
setMethod(f="RunNormalization",
  signature="MultiAssayExperiment",
  definition <- function(object, data, NormMethod){

      groups <- object@metadata$design@List.Factors[object@metadata$design@Factors.Type == "Bio"] %>%
                as.data.frame() %>% unite(col="groups", sep="_")

      coefNorm  = switch(NormMethod,
                        "TMM"=TMM.Normalization(assay(object[[data]]), groups$groups)
                        )
      object@ExperimentList[[data]]@metadata[["Normalization"]] <- list(methode = NormMethod, coefNorm = coefNorm)
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
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param axis The number of PCA axis
#'
#' @exportMethod mvQCdesign
#'
#' @rdname mvQCdesign


setMethod(f="mvQCdesign",
          signature="MultiAssayExperiment",
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
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param axis The number of PCA axis
#'
#' @exportMethod mvQCdata
#' @rdname mvQCdata

setMethod(f="mvQCdata",
          signature="MultiAssayExperiment",
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



#' @title abundanceBoxplot
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param dataType omic data type
#' @exportMethod abundanceBoxplot
#' @rdname abundanceBoxplot
#' 
setMethod(f= "abundanceBoxplot",
          signature = "MultiAssayExperiment",
          definition <- function(object, dataType){

            # this function generate boxplot (abandance distribution) from raw data and normalized data

            #col <- colorPlot(object@design, object@colData, condition="samples")
            sample_names <- row.names(object@colData)
            
            groups  <- object@metadata$design@List.Factors[object@metadata$design@Factors.Type == "Bio"] %>% as.data.frame() %>%
                       unite(col="groups", sep="_", remove = FALSE) %>% mutate(samples=sample_names)

            # normalized data
            pseudo  <- log2(scale(assay(object[[dataType]]), center=FALSE, 
                                  scale=object[[dataType]]@metadata$Normalization$coefNorm$norm.factors)+1) %>% data.table::melt()
    print(object[[dataType]]@metadata$Normalization$coefNorm)
            colnames(pseudo) <- c("feature", "samples", "value")
            pseudo_bis <- full_join(pseudo, groups, by="samples")

            pseudo_bis$samples <- factor(pseudo_bis$samples, levels = unique(pseudo_bis$samples))

            # boxplot
            ggplot(pseudo_bis, aes(x=samples, y=value)) + geom_boxplot(aes(fill=groups)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
              #scale_fill_manual(values=col)
          }
)


#' @title plotPCAnorm
#'
#' @param MultiAssayExperiment
#' @param condition
#' @param color color palette
#'
#' @return
#' @exportMethod plotPCAnorm
#'
#' @examples
setMethod(f= "plotPCAnorm",
          signature = "MultiAssayExperiment",
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
#' @param MultiAssayExperiment
#' @param condition
#' @param colors color palette
#'
#' @return
#' @exportMethod barplotPCAnorm
#' 
#'
#' @examples
setMethod(f= "barplotPCAnorm",
          signature = "MultiAssayExperiment",
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
#' @param MultiAssayExperiment
#' @param threshold
#'
#' @return MultiAssayExperiment
#' @exportMethod FilterLowAbundance
#'
#' @examples
setMethod(f= "FilterLowAbundance",
          signature = "MultiAssayExperiment",
          definition <- function(object, data, threshold){

            objectFilt <- object[[data]]
            
            feature_0  <- objectFilt[rowSums(assay(objectFilt)) <= threshold, ]@NAMES
            
            objectFilt <- objectFilt[rowSums(assay(objectFilt)) > threshold, ]
            
            objectFilt@metadata[["FilteredFeature"]] <-  feature_0
            
            object@ExperimentList[[paste0(data, ".filtred")]] <- objectFilt
            
            return(object)
          })


