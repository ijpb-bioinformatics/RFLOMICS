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



#' @title SetModelMatrixAndContrasts
#' Adapted Function from DicoExpress
#' @param ExpDesign
#'
#' @return
#' @exportMethod SetModelMatrixAndContrasts
#'
#' @examples
#' @examples
#'
setMethod(f="SetModelMatrixAndContrasts",
          signature="ExpDesign",
          definition <- function(object){

  #  => list of all the biological factor

  Factors.Name <- names(object@List.Factors)

  BioFactors <- Factors.Name[which(object@Factors.Type == "Bio")]
  BatchFactors <- Factors.Name[which(object@Factors.Type == "batch")]

  NbBioFactors <- length(BioFactors)
  NbBatchFactors <- length(BatchFactors)

  # => get the model matrix from the formula
  glm_model <- model.matrix(as.formula(object@Model.formula[[1]]),data=as.data.frame( object@List.Factors))
  object@Model.matrix <- glm_model

  #  => list of biological and batch factor from model
  Fact.list <- strsplit(gsub("[ ~]", "",object@Model.formula), split = "+", fixed = TRUE)[[1]]

  BioFactors.formula <- BioFactors[which(BioFactors %in% Fact.list)]
  BatchFactors.formula <- BatchFactors[which(BatchFactors %in% Fact.list)]

  # => Is there any interaction ?

  if(length(grep(":",object@Model.formula))==1){
    Interaction <- TRUE
  }else{
    Interaction <- FALSE
  }

  NbBioFactors.formula <- length(BioFactors.formula)
  NbBatchFactors.formula <- length(BatchFactors.formula)

  ## Contrasts
  FacBio <- 1:(NbBioFactors.formula)
  coeff.name = colnames(glm_model)
  nl = unlist(lapply(object@List.Factors[BioFactors.formula], function(x)
    length(levels(x))))
  nl1 <- levels(object@List.Factors[[BioFactors.formula[1]]])
  contrast.names <- as.character("")

  ## Number of biological factor == 1
  if (NbBioFactors.formula == 1) {
    if ( Interaction == TRUE){
      cat("\n#############################################################\n We can not write Interactions with a single biological factor\n#############################################################\n")
    }
    contrast.factor1 <- matrix(0, ncol = length(coeff.name), nrow = 1)
    for (h in 1:(length(nl1) - 1))
    {
      for (i in (h + 1):length(nl1))
      {
        f2 <- paste(BioFactors.formula[1], nl1[i], sep = "")
        contrast.definition <- rep(0, length(coeff.name))
        f1 <- paste(BioFactors.formula[1], nl1[h], sep = "")
        if (h != 1)
        {
          contrast.definition[which(f1 == coeff.name)] = 1
        }
        contrast.definition[which(f2 == coeff.name)] = (-1)
        contrast.factor1 = rbind(contrast.factor1, contrast.definition)
        contrast.names = c(contrast.names, paste0("[", nl1[h], "-", nl1[i], "]"))
      }
    }
    contrast.matrix <- contrast.factor1

    nbcontrast.type <- dim(contrast.factor1)[1]
    contrast.explanation <- c("effect of biological factor 1")

    colnames(contrast.matrix) = coeff.name
    rownames(contrast.matrix) = contrast.names
  }

  ## Number of biological factor == 2
  if (NbBioFactors.formula == 2) {
    if (Interaction == FALSE) {
      ## Factor 1
      contrast.factor1 <- matrix(0, ncol = length(coeff.name), nrow = 1)
      for (h in 1:(length(nl1) - 1))
      {
        for (i in (h + 1):length(nl1))
        {
          f2 <- paste(BioFactors.formula[1], nl1[i], sep = "")
          contrast.definition <- rep(0, length(coeff.name))
          f1 <- paste(BioFactors.formula[1], nl1[h], sep = "")
          if (h != 1)
          {
            contrast.definition[which(f1 == coeff.name)] = 1
          }
          contrast.definition[which(f2 == coeff.name)] = (-1)
          contrast.factor1 = rbind(contrast.factor1, contrast.definition)
          contrast.names = c(contrast.names, paste0("[", nl1[h], "-", nl1[i], "]"))
        }
      }
      ## Factor 2
      nl2 <- levels(object@List.Factors[BioFactors.formula[2]][[1]])
      contrast.factor2 <- matrix(0, ncol = length(coeff.name), nrow = 1)
      for (j in 1:(length(nl2) - 1))
      {
        for (k in (j + 1):length(nl2))
        {
          f2 <- paste(BioFactors.formula[2], nl2[k], sep = "")
          contrast.definition <- rep(0, length(coeff.name))
          f1 <- paste(BioFactors.formula[2], nl2[j], sep = "")
          if (j != 1)
          {
            contrast.definition[which(f1 == coeff.name)] = 1
          }
          contrast.definition[which(f2 == coeff.name)] = (-1)
          contrast.factor2 = rbind(contrast.factor2, contrast.definition)
          contrast.names = c(contrast.names, paste0("[", nl2[j], "-", nl2[k], "]"))
        }
      }
      contrast.matrix <- unique(rbind(contrast.factor1, contrast.factor2))

      nbcontrast.type <- c(dim(contrast.factor1)[1],dim(contrast.factor2)[1])
      contrast.explanation <- c("effect of biological factor 1","effect of biological factor 2")

      colnames(contrast.matrix) = coeff.name
      rownames(contrast.matrix) = contrast.names
    }

    if ( Interaction == TRUE){
      nl2 <- levels(object@List.Factors[BioFactors.formula[2]][[1]])
      ## Interaction effect between the two biological factors
      contrast.Interaction<-matrix(0,ncol=length(coeff.name),nrow=1)
      for (h in 1:(length(nl1)-1))
      {
        f1<-paste(BioFactors.formula[1],nl1[h],sep="")
        for (i in (h+1):length(nl1))
        {
          f2<-paste(BioFactors.formula[1],nl1[i],sep="")

          for(j in 1:(length(nl2)-1))
          {
            g1<-paste(BioFactors.formula[2],nl2[j],sep="")
            for (k in (j+1):length(nl2))
            {
              g2<-paste(BioFactors.formula[2],nl2[k],sep="")
              contrast.definition<-rep(0,length(coeff.name))
              contrast.definition[grep(paste(f1,g1,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f1,g1,sep=":"),coeff.name),-1,0)
              contrast.definition[grep(paste(f1,g2,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f1,g2,sep=":"),coeff.name),1,0)
              contrast.definition[grep(paste(f2,g1,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f2,g1,sep=":"),coeff.name),1,0)
              contrast.definition[grep(paste(f2,g2,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f2,g2,sep=":"),coeff.name),-1,0)
              contrast.Interaction=rbind(contrast.Interaction,contrast.definition)
              contrast.names=c(contrast.names,paste0("[",nl1[h],"_",nl2[j],"-",nl1[h],"_",nl2[k],"]-[",nl1[i],"_",nl2[j],"-",nl1[i],"_",nl2[k],"]"))
            }
          }
        }
      }
      colnames(contrast.Interaction)=coeff.name
      rownames(contrast.Interaction)=contrast.names

      ## effect of biological factor 1 averaged on biological factor 2
      contrast.factor1<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (h in 1:(length(nl1)-1))
      {
        for (i in (h+1):length(nl1))
        {
          f2<-paste(BioFactors.formula[1],nl1[i],sep="")
          contrast.definition<-rep(0,length(coeff.name))
          f1<-paste(BioFactors.formula[1],nl1[h],sep="")
          if(h!=1)
          {
            contrast.definition[which(f1==coeff.name)]=1
            contrast.definition[grep(paste0(f1,":"),coeff.name)]=(1/nl[2])
          }
          contrast.definition[which(f2==coeff.name)]=(-1)
          contrast.definition[grep(paste0(f2,":"),coeff.name)]=(-1/nl[2])
          contrast.factor1=rbind(contrast.factor1,contrast.definition)
          contrast.names=c(contrast.names,paste0("[",nl1[h],"-",nl1[i],"]"))
        }
      }

      colnames(contrast.factor1)=coeff.name
      rownames(contrast.factor1)=contrast.names


      ## effect of biological factor 2 averaged on biological factor 1
      contrast.factor2<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (j in 1:(nl[2]-1))
      {
        for (k in (j+1):nl[2])
        {
          contrast.definition<-rep(0,length(coeff.name))
          g1<-paste(BioFactors.formula[2],nl2[j],sep="")
          if(j!=1)
          {
            contrast.definition[which(g1==coeff.name)]=1
            contrast.definition[grep(paste0(":",g1),coeff.name)]=(1/nl[1])
          }
          g2<-paste(BioFactors.formula[2],nl2[k],sep="")
          contrast.definition[which(g2==coeff.name)]=(-1)
          contrast.definition[grep(paste0(":",g2),coeff.name)]=(-1/nl[1])
          contrast.factor2=rbind(contrast.factor2,contrast.definition)
          contrast.names=c(contrast.names,paste0("[",nl2[j],"-",nl2[k],"]"))
        }
      }

      colnames(contrast.factor2)=coeff.name
      rownames(contrast.factor2)=contrast.names


      ## effect of biological factor 2 given one level of biological factor 1
      contrast.factor2sachant1<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (h in 1:nl[1])
      {
        h1<-paste(BioFactors.formula[1],nl1[h],sep="")
        for (j in 1:(nl[2]-1))
        {
          for (k in (j+1):nl[2])
          {
            contrast.definition<-rep(0,length(coeff.name))
            g1<-paste(BioFactors.formula[2],nl2[j],sep="")
            if(j!=1)
            {
              contrast.definition[which(g1==coeff.name)]=1
              if(h!=1)
                contrast.definition[grep(paste0(h1,":",g1),coeff.name)]=1
            }
            g2<-paste(BioFactors.formula[2],nl2[k],sep="")
            contrast.definition[which(g2==coeff.name)]=-1
            if(h!=1)
              contrast.definition[grep(paste0(h1,":",g2),coeff.name)]=-1
            contrast.factor2sachant1=rbind(contrast.factor2sachant1,contrast.definition)
            contrast.names=c(contrast.names,paste0("[",nl1[h],"_",nl2[j],"-",nl1[h],"_",nl2[k],"]"))
          }
        }
      }
      colnames(contrast.factor2sachant1)=coeff.name
      rownames(contrast.factor2sachant1)=contrast.names

      ## effect of biological factor 1 given one level of biological factor 2
      contrast.factor1sachant2<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (h in 1:nl[2])
      {
        h1<-paste(BioFactors.formula[2],nl2[h],sep="")
        for (j in 1:(nl[1]-1))
        {
          for (k in (j+1):nl[1])
          {
            contrast.definition<-rep(0,length(coeff.name))
            g1<-paste(BioFactors.formula[1],nl1[j],sep="")
            if(j!=1)
            {
              contrast.definition[which(g1==coeff.name)]=1
              if(h!=1)
                contrast.definition[grep(paste0(g1,":",h1),coeff.name)]=1
            }
            g2<-paste(BioFactors.formula[1],nl1[k],sep="")
            contrast.definition[which(g2==coeff.name)]=(-1)
            if(h!=1)
              contrast.definition[grep(paste0(g2,":",h1),coeff.name)]=(-1)
            contrast.factor1sachant2=rbind(contrast.factor1sachant2,contrast.definition)
            contrast.names=c(contrast.names,paste0("[",nl2[h],"_",nl1[j],"-",nl2[h],"_",nl1[k],"]"))
          }
        }
      }
      colnames(contrast.factor1sachant2)=coeff.name
      rownames(contrast.factor1sachant2)=contrast.names

      contrast.matrix <- unique(rbind(contrast.factor1, contrast.factor2, contrast.factor1sachant2, contrast.factor2sachant1, contrast.Interaction ))
      nbcontrast.type <- unlist(lapply(list(contrast.factor1, contrast.factor2, contrast.factor1sachant2, contrast.factor2sachant1, contrast.Interaction ),
                                       function(x){dim(x)[1]}))
      contrast.explanation <- c("effect of biological factor 1 averaged on biological factor 2" , "effect of biological factor 2 averaged on biological factor 1",
                                "effect of biological factor 2 given one level of biological factor 1", "effect of biological factor 1 given one level of biological factor 2",
                                "Interaction effect between the two biological factors")
    }
  }
  if(nrow(contrast.matrix)==2){
    contrast.matrix <- t(as.data.frame(contrast.matrix[-1,]))
    rownames(contrast.matrix) <- contrast.names[2]
  } else {contrast.matrix <- contrast.matrix[-1,]}


  Contrasts_Names <- data.frame(row.names(contrast.matrix))
  colnames(Contrasts_Names) <- c("Contrasts_Names")


  Contrasts.List <- data.frame("Hypothesis"= row.names(contrast.matrix),
                               "idContrast"=paste("C",1:dim(contrast.matrix)[1],sep=""),
                               "factors"=rep(contrast.explanation,times=nbcontrast.type-1))

  Contrasts.Coeff <-  as.data.frame(t(contrast.matrix))
  names(Contrasts.Coeff) <- Contrasts.List$idContrast

  object@Contrasts.List <- Contrasts.List
  object@Contrasts.Coeff <- Contrasts.Coeff

  return(object)
})

