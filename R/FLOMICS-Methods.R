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
           Model.matrix=vector(),
           Contrasts.List=data.frame(),
           Contrasts.Coeff=data.frame(),
           Contrasts.Sel=vector())
  }

#' @title RunDiffAnalysis
#' @param An object of class [\code{\link{MultiAssayExperiment}]
#' @param data omic data type
#' @param DiffMethod Differential analysis method
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @return MultiAssayExperiment
#' @exportMethod RunDiffAnalysis
#' @examples
#'
setMethod(f="RunDiffAnalysis",
          signature="MultiAssayExperiment",
          definition <- function(object, data, FDR, DiffAnalysisMethod, clustermq){

            # Run the Diff analysis and get the results as a list of object depending of the

            ListOfDiffResults <- switch(DiffAnalysisMethod,
                                     "edgeRglmfit"=edgeR.AnaDiff(object, data, FDR, clustermq)
                                )

            # The results of the methods have to be formated

            # Set an AnaDiff object to
            object@ExperimentList[[data]]@metadata[["AnaDiff"]] <- ListOfDiffResults

            #
            object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]] <- lapply(ListOfDiffResults, function(x){
              res<-topTags(x,10000)
              #DEGs<-res$table[res$table$FDR<=FDR,]
              DEGs<-res$table
              return(DEGs)
            })
            names(object@ExperimentList[[data]]@metadata[["AnaDiffDeg"]]) <- names(ListOfDiffResults)

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
#' @param data data name
#' @param axis The number of PCA axis
#' @param PCA pca axis to plot
#' @param pngFile plot file name to save

#'
#' @exportMethod mvQCdesign
#'
#' @rdname mvQCdesign


setMethod(f="mvQCdesign",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA=c("raw","norm"), axis=5, pngFile=NULL){

            resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <-  object@metadata$colDataStruc["n_dFac"]

            bigdf <- list()
            for(i in 1:n_dFac){
              Factor <-  object@colData[,i]
              nameFac <- names(object@colData)[i]
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
              bigdf[[i]] <- dplyr::bind_rows(df)
            }
            big <- dplyr::bind_rows(bigdf)
            out <- by(data = big, INDICES = big$dfac, FUN = function(m) {
              m <- droplevels(m)
              names(m) <- c("y", "col", "Axis", "dfac", unique(m$dfac), "x")
              m <- ggplot(m,aes(y=y,x=x),aes_string(color=unique(m$dfac)))+
                geom_bar(stat = "identity",position = position_dodge(),aes_string(fill=unique(m$dfac)))+
                facet_grid(as.factor(dfac)~Axis) +
                labs(x = "Samples", y="Coordinates on \n the PCA axis")+
                theme(axis.title.y = element_text(size = 5),
                      axis.title.x=element_text(size = 5),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())
            })
            p <- do.call(grid.arrange, out)
            print(p)
            
            if(! is.null(pngFile)){
              ggsave(filename = pngFile,  plot = p)
            }
})


#' @title multivariate QC data
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param data data name
#' @param axis The number of PCA axis
#' @param PCA pca axis to plot
#' @param pngFile plot file name to save
#'
#' @exportMethod mvQCdata
#' @rdname mvQCdata

setMethod(f="mvQCdata",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA=c("raw","norm"),axis=3, pngFile=NULL){

            resPCA <- object[[data]]@metadata[["PCAlist"]][[PCA]]
            cc <- c(RColorBrewer::brewer.pal(9, "Set1"))

            n_dFac <- object@metadata$colDataStruc["n_dFac"]
            n_qcFac <- dim(object[[data]]@colData[,-c(1:2)])[2]
            var <- names(object[[data]]@colData[,-c(1:2)])

            corA=list()
            VarAxis = list()
            for(i in 1:axis){
            corA[[i]] <- cor(resPCA$ind$coord[,i],
                             as.data.frame(object[[data]]@colData[,-c(1:2)]),
                             method="spearman")
            VarAxis[[i]] <- paste("\n(Var=",round(resPCA$eig[i,2],1),")",sep="")
            }

            df = data.frame("QCparam"=rep(var,axis),
                            "Spearman"= unlist(corA),
                            "Axis"=rep(paste(rep("Axis",axis),1:axis,VarAxis,sep=""), each=n_qcFac))

            p <- ggplot(df,aes(x=Axis, y=abs(Spearman),fill=QCparam))+
              geom_bar(stat="identity",position=position_dodge(),width=0.7)+ylim(0,1)+
              labs(x = "Axis number", y="Cor(Coord_dFactor_PCA,QCparam)")
            
            print(p)
            
            if(! is.null(pngFile)){
              ggsave(filename = pngFile, plot = p)
            }
            
          })



#' @title abundanceBoxplot
#' @param MultiAssayExperiment an object of Class MultiAssayExperiment
#' @param dataType omic data type
#' @param pngFile
#' @exportMethod abundanceBoxplot
#' @rdname abundanceBoxplot
#'
setMethod(f= "abundanceBoxplot",
          signature = "MultiAssayExperiment",
          definition <- function(object, dataType, pngFile=NULL){

            # this function generate boxplot (abandance distribution) from raw data and normalized data

            #col <- colorPlot(object@design, object@colData, condition="samples")
            sample_names <- row.names(object@colData)

            groups  <- object@metadata$design@List.Factors[object@metadata$design@Factors.Type == "Bio"] %>% as.data.frame() %>%
                       unite(col="groups", sep="_", remove = FALSE) %>% mutate(samples=sample_names)

            # normalized data
            pseudo  <- log2(scale(assay(object[[dataType]]), center=FALSE,
                                  scale=object[[dataType]]@metadata$Normalization$coefNorm$norm.factors)+1) %>% reshape2::melt()
            colnames(pseudo) <- c("feature", "samples", "value")
            pseudo_bis <- full_join(pseudo, groups, by="samples")

            pseudo_bis$samples <- factor(pseudo_bis$samples, levels = unique(pseudo_bis$samples))

            # boxplot
            p <- ggplot(pseudo_bis, aes(x=samples, y=value)) + geom_boxplot(aes(fill=groups)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
              xlab(paste0(dataType, " samples"))
              
              #scale_fill_manual(values=col)
            print(p)
            
            if(! is.null(pngFile)){
              ggsave(filename = pngFile, plot = p)
            }
            
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
          definition <- function(object, data, PCA, PCs=c(1,2), condition="groups", pngFile){

            #
            PC1 <- paste("Dim.",PCs[1], sep="")
            PC2 <- paste("Dim.",PCs[2], sep="")

            #col <- colorPlot(object@design, object@colData, condition=condition)

            #sample_names <- row.names(object@colData)

            groups    <- object@metadata$design@List.Factors[object@metadata$design@Factors.Type == "Bio"] %>% as.data.frame() %>%
                         unite(col="groups", sep="_", remove = FALSE) #%>% mutate(samples=sample_names)


            factors   <- object@metadata$design@List.Factors %>% as.data.frame() %>%
                         unite(., col="samples", sep="_", remove = FALSE) %>% mutate(groups = groups$groups)

            score     <- FlomicsMultiAssay[[data]]@metadata$PCAlist[[PCA]]$ind$coord[, PCs] %>% as.data.frame() %>%
                         mutate(samples=row.names(.)) %>% full_join(., factors, by="samples")

            var1 <- round(FlomicsMultiAssay[[data]]@metadata$PCAlist[[PCA]]$eig[PCs,2][1], digit=3)
            var2 <- round(FlomicsMultiAssay[[data]]@metadata$PCAlist[[PCA]]$eig[PCs,2][2], digit=3)

            
            p <- ggplot(score, aes_string(x=PC1, y=PC2, color=condition))  +
              geom_point(size=3) +
              geom_text(aes(label=samples), size=3, vjust = 0) +
              xlab(paste(PC1, " (",var1,"%)", sep="")) +
              ylab(paste(PC2, " (",var2,"%)", sep="")) +
              geom_hline(yintercept=0, linetype="dashed", color = "red") +
              geom_vline(xintercept=0, linetype="dashed", color = "red") +
              theme(strip.text.x = element_text(size=8, face="bold.italic"),
                    strip.text.y = element_text(size=8, face="bold.italic")) #+
              #scale_color_manual(values=col$colors)
            
            print(p)
            ggsave(filename = pngFile,  plot = p)
            
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



#' @title RunPCA
#' @param MultiAssayExperiment
#' @param data omic data type
#' @param PCA norm or raw
#' @return MultiAssayExperiment
#' @exportMethod RunPCA
#' @examples
#'
setMethod(f="RunPCA",
          signature="MultiAssayExperiment",
          definition <- function(object, data, PCA){

            if(PCA=="raw"){
            pseudo <- log2(scale(assay(object[[data]]),
                                      center=FALSE)+1)
            }
            else if(PCA=="norm"){
            pseudo <- log2(scale(assay(object[[data]]),
                                        center=FALSE,
                                        scale=object[[data]]@metadata[["Normalization"]]$coefNorm$norm.factors)+1)
            }

            pca <- FactoMineR::PCA(t(pseudo),ncp = 5,graph=F)

            object@ExperimentList[[data]]@metadata[["PCAlist"]][[PCA]]<- pca

            return(object)
          }
)



#' @title SetModelMatrix
#'
#' @param ExpDesign 
#'
#' @return
#' @exportMethod SetModelMatrix
#'
#' @examples
setMethod(f="SetModelMatrix",
          signature="ExpDesign",
          definition <- function(object){
            
            #  => list of fact from model
            Fact.list <- strsplit(gsub("[ ~]", "",object@Model.formula), split = "+", fixed = TRUE)[[1]]
            Fact.list <- Fact.list[str_detect(Fact.list, pattern=":", negate = TRUE)]
            
            #  => list of fact bio
            Fact.Bio  <- names(object@Factors.Type[which(object@Factors.Type == "Bio")])
            
            #  => factor vectors to set model.matrix `
            #  * stock ref values
            ref.list <- vector()
            for (i in Fact.list) {
              assign(i, object@List.Factors[[i]])
              ref.list <- c(ref.list,levels(object@List.Factors[[i]])[1])
            }
            
            #  => set model matrix
            #  * change colnames
            model.design.matrix <- stats::model.matrix(as.formula(object@Model.formula))
            design.colnames <- colnames(model.design.matrix)
            #  colnames
            for (i in names(object@List.Factors)) {
              design.colnames <- gsub(i, "", design.colnames)
            }
            design.colnames <- gsub(":", "_", design.colnames)
            colnames(model.design.matrix) <- design.colnames
            object@Model.matrix <- model.design.matrix

            # => Get all the contrasts
            contrasts <- list()
            contrast.coef <- list()
            contrasts.nbr <- 0
            for(i in Fact.list[which(Fact.list %in% Fact.Bio)]){
              
              mat <- as.vector(unique(object@List.Factors[[i]])) %>% utils::combn(.,2)
              
              contrast.tmp <- list()
              
              contrast.tmp[["hypoth"]] <- apply (mat, 2, function(x) {
                paste(x, collapse=" - ")
              })
              
              contrast.tmp[["hypoth_tmp"]] <- contrast.tmp[["hypoth"]]
              
              for (ref in ref.list) {
                contrast.tmp[["hypoth_tmp"]] <- gsub(paste(ref, "-"), "", contrast.tmp[["hypoth_tmp"]])
                contrast.tmp[["hypoth_tmp"]] <- gsub(paste("-", ref), "", contrast.tmp[["hypoth_tmp"]])
              }
              
              contrast.tmp[["idContrast"]] <- paste("C", (contrasts.nbr+1):(contrasts.nbr+length(contrast.tmp[["hypoth"]])), sep = "")
              
              contrasts.nbr <- contrasts.nbr+length(contrast.tmp[["hypoth"]])
              contrasts[[i]] <- data.frame(contrast.tmp) %>% mutate(factors=i)
                   
              ## contrast coef
              #contrast.coef[[i]] <- makeContrasts(contrasts = contrasts[[i]]$hypoth_tmp, levels = model.design.matrix) %>% as.data.frame()
              #colnames(contrast.coef[[i]]) <- contrasts[[i]]$idContrast
            }
            object@Contrasts.List  <- contrasts %>% purrr::reduce(rbind)
            
            ## contrast coef
            object@Contrasts.Coeff <- makeContrasts(contrasts = object@Contrasts.List$hypoth_tmp, levels = model.design.matrix) %>% as.data.frame()
            colnames(object@Contrasts.Coeff) <- object@Contrasts.List$idContrast
            
            return(object)

})

