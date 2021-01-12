#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' GetDesignFromNames
#'
#' @param samples_name a vector of sample names giving the designs factor, each
#' separated by "_"
#'
#' @return a named list of two elements
#' * nb_dfac: the number of design factors
#' * tmpDesign: a data frame with the names of design factors in columns (ie: dFac1) and all
#' the modalities in row
#' @md
#' @export
#'
#' @examples
#'
#'
GetDesignFromNames <- function(samples_name){

  # Get the number of design factor and the factors from the names of the count matrix
  nb_dFac <- stringr::str_count(samples_name,pattern="_")+1
  # Test if the number of factor are the same for all sample names
  try(if(var(nb_dFac) != 0 ) stop("Column names do not have the same level of factor"))
  nb_dFac <- nb_dFac[1]
  names(nb_dFac) <- "n_dFac"

  #  Get all the factors from the names of the count matrix
  tmpDesign <- tibble::tibble(design=samples_name) %>%
    tidyr::separate(.,design,into=paste("dFac",1:nb_dFac["n_dFac"],sep="")) %>%
    dplyr::mutate_all(.,as.factor)

  return(list("nb_dFac"=nb_dFac,"tmpDesign"=tmpDesign))
}

#' GetModelFormulae
#'
#' @param Factors.Name
#' @param Factors.Type
#'
#' @return a named list of formulae
#' @export
#'
#' @examples
#'
GetModelFormulae <- function(Factors.Name,Factors.Type){


  formulae <- list()

  FacBio <- Factors.Name[which(Factors.Type == "Bio")]
  FacBatch <- Factors.Name[which(Factors.Type == "batch")]

  nFac <- length(FacBio)

  getF <- function(x,FacBatch){
    update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+","(",paste(x,collapse="+"),")^2")),new=~.)
  }
  getF2 <- function(x,FacBatch){
    update(as.formula(paste(paste("~ ",FacBatch,collapse ="+"),"+",paste(x,collapse="+"))),new=~.)
  }

  for(i in 1:nFac){

    formulae <- c(formulae, apply(combn(FacBio,i),2,getF, FacBatch=FacBatch))
    if(i !=1){
      formulae <- c(formulae, apply(combn(FacBio,i),2,getF2, FacBatch=FacBatch))
    }
    #formulae[[i]]  <- apply(combn(FacBio,i),2,getF, FacBatch=FacBatch)
  }

  #formulae[[nFac+1]]<-apply(combn(FacBio,i),2,getF2,FacBatch=FacBatch)

  formulae <- unlist(formulae)
  names(formulae) <- unlist(as.character(formulae))

  #formulae[[nFac+1]]<- apply(combn(FacBio,i),2,getF2,FacBatch=FacBatch)

  #formulae <- unlist(formulae)
  #names(formulae) <- unlist(as.character(formulae))

  return(formulae)
}


#' GetContrasts
#'
#' @param An object of class design
#'
#' @return An object of class design
#' @export
#'
#' @examples
#'
GetContrasts <- function(design){
  # function qui est en cours d'ecriture par Christine
  data.frame()
}

#' @title TMM.Normalization
#' Interface to the calcNormFactors functionof the edgeR package  with the choosen TMM parameters as the normalization method
#' @param counts numeric matrix of read counts
#' @param groups vector or factor giving the experimental group/condition for each sample/library.
#' @return a data.frame with a row for each sample and columns group, lib.size and norm.factors containing the group labels, library sizes and normalization factors. Other columns can be optionally added to give more detailed sample information.
#' @export
#'
#' @examples
TMM.Normalization <- function(counts, groups){
  dge <- edgeR::DGEList(counts=counts, group=groups)
  dge <- edgeR::calcNormFactors(dge,method="TMM")
  nf  <- dge$samples
  return(nf)
}


#' @title edgeR.AnaDiff
#'
#' @param object an object of class [\code{\link{MultiAssayExperiment}]
#' @param data Omic data type
#' @param clustermq A boolean indicating if the constrasts have to be computed in local or in a distant machine
#' @return A list of object of class [\code{\link{DGELRT}]
#' @export
#'
#' @examples
edgeR.AnaDiff <- function(object, data, clustermq){
  
  # retrieve the design matrix
  model_matrix <- model.matrix(as.formula(object@metadata$design@Model.formula),
                               data=as.data.frame(object@metadata$design@List.Factors))

  # Construct the DGE obect
  dge <- edgeR::DGEList(counts=assay(object@ExperimentList[[data]]),
                        group=object@ExperimentList[[data]]@metadata$Normalization$coefNorm$group,
                        lib.size =object@ExperimentList[[data]]@metadata$Normalization$coefNorm$lib.size,
                        norm.factors = object@ExperimentList[[data]]@metadata$Normalization$coefNorm$norm.factors)

  # Run the model
  print("[cmd] dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)")
  dge <- edgeR::estimateGLMCommonDisp(dge, design=model_matrix)
  print("[cmd] dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)")
  dge <- edgeR::estimateGLMTrendedDisp(dge, design=model_matrix)
  print("[cmd] dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)")
  dge <- edgeR::estimateGLMTagwiseDisp(dge, design=model_matrix)
  print("[cmd] fit.f <- edgeR::glmFit(dge,design=model_matrix)")
  fit.f <- edgeR::glmFit(dge,design=model_matrix)

  # test clustermq

  if(clustermq == TRUE){

     # Fonction to run on contrast per job
     # y is the model, Contrasts are stored in a matrix, by columns
    
      fx <- function(x){
        edgeR::glmLRT(y, contrast = unlist(z[x,]))
      }
    
      ListRes <- clustermq::Q(fx, x=1:length(object@metadata$design@Contrasts.Sel$contrast),
        export=list(y=fit.f,z=object@metadata$design@Contrasts.Coeff),
        n_jobs=length(object@metadata$design@Contrasts.Sel$contrast),pkgs="edgeR")
      
  }
  else{
    
     ListRes <-  lapply(object@metadata$design@Contrasts.Sel$contrast, function(x){

       edgeR::glmLRT(fit.f, contrast = unlist(object@metadata$design@Contrasts.Coeff[x,]))
       
     })
  }

  names(ListRes) <- object@metadata$design@Contrasts.Sel$contrastName
  return(ListRes)
}



#' @title colorPlot
#'
#' @param design
#' @param ColData
#' @param condition
#'
#' @return
#' @export
#'
#' @examples
colorPlot <- function(design, ColData, condition="samples"){

  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

    if(condition == "samples"){

      # combine only bio fact
      groups <- design@List.Factors[design@Factors.Type == "Bio"] %>%
        as.data.frame() %>% unite(col="groups", sep="_")
      list.cond <- factor(groups$groups)
    }
    else{

      list.cond <- design@List.Factors[[condition]]
    }

    len.cond  <- length(levels(list.cond))

    colors    <- getPalette(len.cond)

    #col <- colors[list.cond]
    col <- colors[levels(list.cond)]
    #names(col) <- row.names(ColData)
    names(col) <- row.names(levels(list.cond))

    return(col)
}


#' plotLibSize
#'
#' @param abundances
#' @param dataName
#' @param pngFile 
#' @return plot
#' @export
#'
#' @examples
plotLibSize <- function(abundances, dataName, pngFile=NULL){

  samples     <- colnames(abundances)
  libSizeNorm <- data.frame ( value = colSums(abundances, na.rm = TRUE) , samples=samples)

  libSizeNorm$samples <- factor(libSizeNorm$samples, levels = libSizeNorm$samples)

  p <- ggplot(libSizeNorm, aes(x=samples,y=value, fill=samples)) + geom_bar( stat="identity" ) +
    xlab(paste0(dataName, " samples")) + ylab("Library Size") + 
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          legend.position  = "none")
          #axis.text.x     = element_blank(),
          #axis.ticks      = element_blank())
          #legend.key.size = unit(0.3, "cm"))
          #legend.text     = element_text(size=5))
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }
}


#' plotDistr
#'
#' @param abundances matrix or dataframe of feature/gene abundances/counts
#' @param dataName name of dataset
#' @param pngFile png file name
#' @return plot
#' @export
#'
#' @examples plotDistr(assay(MAE), "dataset1", "tmp/countDist.png")
plotDistr <- function(abundances, dataName, pngFile=NULL){

  pseudo_counts <- log2(abundances+1) %>% reshape2::melt()
  colnames(pseudo_counts) <- c("features", "samples", "value")

  p <- ggplot(pseudo_counts) +
    geom_density(aes(value, color=samples) ) +
    xlab(paste0(dataName, " log2(feature abundances)")) + 
    theme(legend.position='none')
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }
}




#' pvalue.plot
#'
#' @param data 
#' @param contrast
#' @param pngFile 
#' @return plot
#' @export
#'
#' @examples
pvalue.plot <- function(data, contrast, pngFile=NULL){

  p <- ggplot(data=data) + geom_histogram(aes(x=PValue), bins = 200)
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }

}







#' Plot the balance of data in an experimental design
#' 
#' This function provides easy visualization of the balance of data in a data set given a specified experimental design. This function is useful for identifying 
#' missing data and other issues. The core of this function is from the function ezDesign in the package ez.
#' 
#' @param counts : the number of data in each cell of the design 
#' @param cell_border_size : Numeric value specifying the size of the border seperating cells (0 specifies no border)
#'
#' @return A printable/modifiable ggplot2 object.
#' @export
#'
#' @author Christine Paysant-Le Roux
plotExperimentalDesign <- function(counts, cell_border_size = 10){
  if (names(counts)[ncol(counts)] != "Count"){
    stop("the last column of the input data frame must be labelled Count")
  }
  if(ncol(counts) == 5){
    x <- names(counts)[1]
    y <- names(counts)[2]
    row <- names(counts)[3]
    col <- names(counts)[4]
  } else if(ncol(counts) == 4){
    x <- names(counts)[1]
    y <- names(counts)[2]
    row <- names(counts)[3]
    col <- NULL
  } else if(ncol(counts) == 3){
    x <- names(counts)[1]
    y <- names(counts)[2]
    row <- NULL
    col <- NULL
  } else if(ncol(counts) == 2){
    x <- NULL
    y <- names(counts)[1]
    row <- NULL
    col <- NULL
  } else {
    stop("data frame with less than 2 columns")
  }
  # rename two first column names of counts with x and y 
  if(!is.null(x)){
    x_lab <- names(counts)[names(counts)==x]
    names(counts)[names(counts)==x] <- 'x'
  } else {
    x_lab <- ""
    counts$x <- rep(1, nrow(counts))
  }
  
  y_lab = names(counts)[names(counts)==y]
  names(counts)[names(counts)==y] <- 'y'
  # get the levels of one factor
  getFactorLevels <- function(counts, factorVar){
    if(!is.numeric(counts[[factorVar]])){
      counts[[factorVar]] = factor(counts[[factorVar]])
      x_vals = as.character(levels(counts[[factorVar]]))
    }else{
      x_vals = as.character(sort(unique(counts[[factorVar]])))
    }
    return(x_vals)
  }
  x_vals <- getFactorLevels(counts, "x")
  # recode factor with integer (1 for the first level, 2 for the second level, etc)
  counts$x = as.numeric(factor(counts$x))
  y_vals <- getFactorLevels(counts, "y")
  counts$y = as.numeric(factor(counts$y))
  # cell border size
  if(length(unique(counts$y))>length(unique(counts$x))){
    cell_border_size = cell_border_size/length(unique(counts$y))
  }else{
    cell_border_size = cell_border_size/length(unique(counts$x))
  }
  # for row or col (for faceting)
  # replace colname with variable name and each level with the paste of the factor name, =, the factor level
  replaceRowOrColIntoCounts <- function(counts, facetingVariable, newColumnName){
    if(!is.factor(counts[,names(counts)== facetingVariable])){
      counts[,names(counts) == facetingVariable] = factor(counts[,names(counts) == facetingVariable])
    }
    levels(counts[,names(counts) == facetingVariable]) = paste(
      names(counts)[names(counts) == facetingVariable]
      , levels(counts[,names(counts) == facetingVariable])
      , sep = ' = '
    )
    names(counts)[names(counts) == facetingVariable] = newColumnName
    return(counts)
  }
  
  if(!is.null(row)){
    counts <- replaceRowOrColIntoCounts (counts, row, "row")
  }
  if(!is.null(col)){
    counts <- replaceRowOrColIntoCounts (counts, col, "col")
  }
  # calculate parameters for geom_rect
  # xmin - (required) left edge of rectangle
  # xmax - (required) right edge of rectangle
  # ymin - (required) bottom edge of rectangle
  # ymax - (required) top edge of rectangle 
  counts$ymin = counts$y-.5
  counts$ymax = counts$y+.5
  counts$xmin = counts$x-.5
  counts$xmax = counts$x+.5
  
  counts$Count <- as.factor(counts$Count)
  
  library(ggplot2)
  p <- ggplot2::ggplot( data = counts, aes_string(ymin = 'ymin', ymax = 'ymax', xmin = 'xmin', xmax = 'xmax' , fill = 'Count')) + geom_rect() + labs(x=x_lab,y=y_lab)
  p <- p + theme(
    panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , legend.background = element_rect(colour='transparent',fill='transparent')
  )
  if(cell_border_size>0){
    p = p + geom_rect(
      size = cell_border_size
      , colour = 'grey90'
      , show.legend = FALSE
    )
  }
  p = p + scale_x_continuous(
    breaks = sort(unique(counts$x))
    , labels = x_vals
  )
  p = p + scale_y_continuous(
    breaks = sort(unique(counts$y))
    , labels = y_vals
  )
  if(!is.null(row)){
    if(!is.null(col)){
      p = p + facet_grid(row~col)
    }
    else{
      p = p + facet_grid(row~.)
    }
  }else{
    if(!is.null(col)){
      p = p + facet_grid(.~col)
    }
  }
  return(p)
}




#' SummarizeExpDesignInput
#'
#' @param ExpDesign data.frame 
#' @param 
#' @param 
#' @return conditSumm data.frame
#' @export
#'
#' @examples
SummarizeExpDesignInput <- function(ExpDesign){
  
  condit_list <- lapply(colnames(ExpDesign[,-1]), function(x){
    conditCount <- table(as.data.frame(ExpDesign)[,x])
    paste(names(conditCount), conditCount, sep=" x ") %>% paste(., collapse =" ; ")
    })
  
  #names(condit_list) <- colnames(ExpDesign[,-1])
  conditSumm <- data.frame("Factors"    = colnames(ExpDesign[,-1]),
                           "Conditions" = condit_list %>% unlist ,
                           "Comments"   = lapply(colnames(ExpDesign[,-1]), function(x){if_else(length(unique(table(as.data.frame(ExpDesign)[,x]))) == 1, "OK", "WARNING") }) %>% unlist)
  
  return(conditSumm)
}




######################################## get contrast ##################################################

# it is possible to define contrast combinations that are specifically suited to a particular experimental design and set of research questions

#Contrasts are used to test whether the levels of an effect are significantly different from one another. You
#can specify a contrast for each factor in the model. Contrasts represent linear combinations of the
#parameters.

# Comparison (or contrast) procedures are used to test more specific hypotheses about differences between means
# comparisons of the two drug groups to the control group (  a complex comparison)
# The contrasts are formed by applying a set of weights, called contrast coefficients, to the means
# A contrast is a test of the difference between the means of two groups from the ANOVA. There are two categories of contrasts among the groups tested by ANOVA, simple and complex.
# A simple contrast is a test of the difference between any two pairs, such as Experimental Group 1 and Control Group 2. A complex contrast is a test of the difference between combinations of groups.
# An example of a complex contrast is a test of the difference between a subgroup created by combining Experimental Groups 1, 2 and 4 combined, and a subgroup created by combining Control Groups 1 and 3

# pairwise comparisons between the groups
# biological_factors_list, interaction_in_model, model_matrix

## functions to define part of simple contrasts
#' define part of simple contrast data frame
#'
#' @param treatmentFactorsList
#' @param i
#' @param j
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
define_partOfSimpleContrast_df <- function (treatmentFactorsList, i, j) {
  comparisonPart <- treatmentFactorsList
  # combn(x,2) generate all combinations of the elements of x taken 2 at a time
  vectorFromCombn <- combn(treatmentFactorsList[[i]],2)[j,]
  comparisonPart[[i]] <- vectorFromCombn
  df_comparisonPart <- expand.grid(comparisonPart)
  data.table::setDT(df_comparisonPart)
  
  # paste all the column of the data table
  #df_comparisonPart[, contrastPart := do.call(paste, c(.SD, sep = "_")), .SDcols = names(df_comparisonPart)]
  df_comparisonPart <- df_comparisonPart %>% tidyr::unite(contrastPart, sep="_", remove=F) 
  
  #colnameFactor_i <- names(df_comparisonPart)[i]
  colnameFactor_i <- names(treatmentFactorsList)[i]
  
  #df_comparisonPart[, comparisonPart := df_comparisonPart[[colnameFactor_i]]]
  df_comparisonPart <- df_comparisonPart %>% dplyr::mutate(comparisonPart = df_comparisonPart[[colnameFactor_i]])
  
  if(length(names(treatmentFactorsList)) != 1){
    colnamesToKeep <- setdiff(names(df_comparisonPart),c("contrastPart", "comparisonPart", colnameFactor_i))
    #df_comparisonPart[, fixFactor := do.call(paste, c(.SD, sep = "_")), .SDcols = colnamesToKeep]
    df_comparisonPart <- df_comparisonPart %>% tidyr::unite(fixFactor, all_of(colnamesToKeep), sep = "_", remove = F)
  }else{
    df_comparisonPart <- df_comparisonPart %>% dplyr::mutate(fixFactor= NA)
  }
  
  
  colnamesToDelete <- names(treatmentFactorsList)
  #df_comparisonPart[, (colnamesToDelete) := NULL]
  df_comparisonPart <- df_comparisonPart %>% dplyr::select(-all_of(colnamesToDelete))
  
  nameColumnContrast <- paste0("contrastPart", j)
  nameColumnComparison <- paste0("comparisonPart", j)
  setnames(df_comparisonPart, c("contrastPart", "comparisonPart"), c(nameColumnContrast, nameColumnComparison))
  return(df_comparisonPart)
}
#' compute a data table with all pairwise comparisons of one factor
#'
#' @param treatmentFactorsList
#' @param i
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
simpleContrastForOneFactor <- function (treatmentFactorsList, i){
  df_FirstComparisonPart <- define_partOfSimpleContrast_df(treatmentFactorsList,i,2)
  #df_FirstComparisonPart[,fixFactor := NULL]
  df_FirstComparisonPart <- df_FirstComparisonPart %>% dplyr::select(-fixFactor)
  
  df_SecondComparisonPart <- define_partOfSimpleContrast_df(treatmentFactorsList,i,1)
  df_simpleContrasts_factor <- cbind(df_FirstComparisonPart, df_SecondComparisonPart)
  
  #df_simpleContrasts_factor[, contrast := paste0("(", contrastPart2, " - ", contrastPart1, ")")]
  #df_simpleContrasts_factor[, groupComparison := paste0("(", comparisonPart2, " - ", comparisonPart1, ")")]
  
  df_simpleContrasts_factor <- df_simpleContrasts_factor %>% 
    dplyr::mutate(contrast        = paste0("(", contrastPart2,   " - ", contrastPart1, ")"),
                  groupComparison = paste0("(", comparisonPart2, " - ", comparisonPart1, ")"))
  
  
  # case where fixFactor column is empty (NA inside)
  emptycolFixFactor <- unique(is.na(df_simpleContrasts_factor$fixFactor))
  if(emptycolFixFactor){
    #df_simpleContrasts_factor[, contrastName := groupComparison]
    df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::mutate(contrastName = groupComparison)
    
  } else {
    #df_simpleContrasts_factor[, contrastName := paste0(groupComparison, " in ", fixFactor )]
    df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::mutate(contrastName = paste0(groupComparison, " in ", fixFactor ))
  }
  #df_simpleContrasts_factor[, type := "simple"]
  df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::mutate(type = "simple")
  
  colnamesToDelete <- c("contrastPart2", "comparisonPart2", "contrastPart1", "comparisonPart1")
  
  #df_simpleContrasts_factor[, (colnamesToDelete) := NULL]
  #setcolorder(df_simpleContrasts_factor, c(names(df_simpleContrasts_factor)[2:length(names(df_simpleContrasts_factor))], names(df_simpleContrasts_factor)[1]))
  df_simpleContrasts_factor <- df_simpleContrasts_factor %>% dplyr::select(-all_of(colnamesToDelete)) %>%
    dplyr::select("contrast", "groupComparison", "contrastName", "type", "fixFactor")
  
  
  
  return(df_simpleContrasts_factor)
}

#' define all simple contrasts
#'
#' @param treatmentFactorsList
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
defineAllSimpleContrasts <- function(treatmentFactorsList){
  # create a data table with 5 columns
  # empty data table
  allSimpleContrast_df <- data.table(contrast = character(), groupComparison = factor(), contrastName = character(), type = character(), fixFactor = factor())
  # create each data frame and rbind it to the allSimpleContrast_df
  for(i in seq_along(treatmentFactorsList)){
    dataTableToCreate <- simpleContrastForOneFactor(treatmentFactorsList, i)
    allSimpleContrast_df <- rbind(allSimpleContrast_df, dataTableToCreate)
  }
  #allSimpleContrast_df[,contrastCoeff := sapply(contrast, function(x) defineCoefficient(x, colnamesGLMdesign))]
  return(allSimpleContrast_df[])
}

#' Define averaged contrasts
#'
#' @param allSimpleContrast_df : a data frame with all the simple contrasts comparisons (test of the difference between any two pairs of groups)
#'
#' @return allAveragedContrasts_df : a data frame with all the averaged contrasts
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
define_averaged_contrasts <- function(allSimpleContrast_df){
  #allAveragedContrasts_df <- allSimpleContrast_df[, list(contrast = paste0(paste0(paste0("(", paste(contrast, collapse=" + ")),")/"),.N), meanIn = paste(fixFactor, collapse=" + ")), by = groupComparison]
  #allAveragedContrasts_df[, type :="mean"]
  allAveragedContrasts_df <- allSimpleContrast_df %>% group_by(groupComparison) %>% dplyr::add_tally() %>% 
    dplyr::mutate(contrast= paste0(paste0("(", paste(contrast, collapse=" + ")),")/", n), 
           meanIn  = paste(fixFactor, collapse=" + "),
           type    = "mean") %>% 
    dplyr::select(-contrastName, -fixFactor, -n) %>% unique() %>% data.table()

  #allAveragedContrasts_df[, contrastName := paste(groupComparison, "mean", sep = " in ")]
  allAveragedContrasts_df <- allAveragedContrasts_df %>% dplyr::mutate(contrastName = paste(groupComparison, "mean", sep = " in "))
  
  setcolorder(allAveragedContrasts_df, c("contrast", "groupComparison", "contrastName", "type", "meanIn"))
  #  allAveragedContrasts_df <- allSimpleContrast_df[, list(meanIn = paste(fixFactor, collapse=" + ")), by = groupComparison]
  return(allAveragedContrasts_df[])
}
# define contrasts for interactions

#' define a data frame with part of interaction contrast
#'
#' @param treatmentFactorsList
#' @param i
#' @param j
#' @param k
#' @param row_i
#' @param row_j
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
define_partOfInteractionContrast_df <- function (treatmentFactorsList, i, j, k, row_i, row_j) {
  comparisonPart <- treatmentFactorsList
  # combn(x,2) generate all combinations of the elements of x taken 2 at a time
  comparisonPart [[i]] <- combn(treatmentFactorsList[[i]],2)[row_i,]
  comparisonPart [[j]] <- combn(treatmentFactorsList[[j]],2)[row_j,]
  df_comparisonPart <- expand.grid(comparisonPart)
  data.table::setDT(df_comparisonPart)
  # paste all the column of the data table
  #df_comparisonPart[, contrastPart := do.call(paste, c(.SD, sep = "_")), .SDcols = names(df_comparisonPart)]
  df_comparisonPart <- df_comparisonPart %>% tidyr::unite(contrastPart_bis, names(df_comparisonPart), sep="_", remove=F) %>%
    dplyr::mutate(contrastPart = contrastPart_bis) %>% dplyr::select(-contrastPart_bis)
  
  colnameFactor_i <- names(df_comparisonPart)[i]
  colnameFactor_j <- names(df_comparisonPart)[j]
  
  #data.table::setDT(df_comparisonPart)
  #df_comparisonPart[, comparisonPart := df_comparisonPart[[colnameFactor_i]]]
  #df_comparisonPart[, fixPart := df_comparisonPart[[colnameFactor_j]]]
  df_comparisonPart <- df_comparisonPart %>% 
    dplyr::mutate(comparisonPart = df_comparisonPart[[colnameFactor_i]]) %>%
    dplyr::mutate(fixPart        = df_comparisonPart[[colnameFactor_j]])
  
  
  colnamesToKeep <- setdiff(names(treatmentFactorsList),c(colnameFactor_i, colnameFactor_j))
  #df_comparisonPart[, outsideGroup := do.call(paste, c(.SD, sep = "_")), .SDcols = colnamesToKeep]
  if(length(colnamesToKeep)){
    
    df_comparisonPart <- df_comparisonPart %>% tidyr::unite(outsideGroup_bis, all_of(colnamesToKeep), sep="_", remove=F) %>% 
      dplyr::mutate(outsideGroup = outsideGroup_bis) %>% dplyr::select(-outsideGroup_bis)
  }else{
    
    df_comparisonPart <- df_comparisonPart %>% dplyr::mutate(outsideGroup = NA)
  }
  
  colnamesToKeep <- setdiff(names(df_comparisonPart),c("contrastPart", "comparisonPart", "fixPart", "outsideGroup", colnameFactor_i))
  #df_comparisonPart[, fixFactor := do.call(paste, c(.SD, sep = "_")), .SDcols = colnamesToKeep]
  df_comparisonPart <- df_comparisonPart %>% tidyr::unite(fixFactor_bis, all_of(colnamesToKeep), sep="_", remove=F) %>% 
    dplyr::mutate(fixFactor = fixFactor_bis) %>% dplyr::select(-fixFactor_bis)
  
  colnamesToDelete <- names(treatmentFactorsList)
  #df_comparisonPart[, (colnamesToDelete) := NULL]
  df_comparisonPart <- df_comparisonPart %>% dplyr::select(-all_of(colnamesToDelete))
  
  
  nameColumnContrast <- paste0("contrastPart", k)
  nameColumnComparison <- paste0("comparisonPart", k)
  nameFixFactor <- paste0("fixFactor", k)
  namePartFixFactor <- paste0("fixPart", k)
  nameOutsideGroup <- paste0("outsideGroup", k)
  setnames(df_comparisonPart, c("contrastPart", "comparisonPart", "fixFactor", "fixPart", "outsideGroup"),
           c(nameColumnContrast, nameColumnComparison, nameFixFactor, namePartFixFactor, nameOutsideGroup))
  return(df_comparisonPart)
}
#' define interaction constrast for pairs of biological factors
#'
#' @param treatmentFactorsList
#' @param i
#' @param j
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
defineInteractionConstrastForPairsOfFactors <- function(treatmentFactorsList, i, j){
  df_part1<- define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 1, 2, 2)
  df_part2<- define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 2, 1, 2)
  df_part3<- define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 3, 2, 1)
  df_part4<- define_partOfInteractionContrast_df (treatmentFactorsList, i, j, 4, 1, 1)
  df_interactionContrasts <- cbind(df_part1, df_part2, df_part3, df_part4)
  #df_interactionContrasts[, contrast := paste0("(", "(", contrastPart1, " - ", contrastPart2, ")"," - ", "(", contrastPart3, " - ", contrastPart4, ")", ")")]
  #df_interactionContrasts[, groupComparison := paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " vs ", "(", fixPart1, " - ", fixPart3, ")")]
  #df_interactionContrasts[, contrastName := paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " in ", fixFactor1, " - ", "(", comparisonPart3, " - ", comparisonPart4, ")", " in ", fixFactor3 )]
  #df_interactionContrasts[, type := "interaction"]
  df_interactionContrasts <- df_interactionContrasts %>% dplyr::mutate(contrast = paste0("(", "(", contrastPart1, " - ", contrastPart2, ")"," - ", 
                                                                                              "(", contrastPart3, " - ", contrastPart4, ")", ")"),
                                                                       groupComparison = paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " vs ", 
                                                                                                "(", fixPart1, " - ", fixPart3, ")"),
                                                                       contrastName  = paste0("(", comparisonPart1, " - ", comparisonPart2, ")", " in ", fixFactor1, " - ", 
                                                                                              "(", comparisonPart3, " - ", comparisonPart4, ")", " in ", fixFactor3 ),
                                                                       type = "interaction")

  colnamesToDelete <- c("contrastPart1",  "comparisonPart1", "fixFactor1", "fixPart1", "outsideGroup1",
                        "contrastPart2", "comparisonPart2", "fixFactor2", "fixPart2", "outsideGroup2",
                        "contrastPart3", "comparisonPart3", "fixFactor3", "fixPart3", "outsideGroup3",
                        "contrastPart4", "comparisonPart4", "fixFactor4", "fixPart4")
  #df_interactionContrasts[, (colnamesToDelete) := NULL]
  df_interactionContrasts <- df_interactionContrasts %>% dplyr::select(-all_of(colnamesToDelete))
    
  setnames(df_interactionContrasts, "outsideGroup4", "outsideGroup")
  #df_interactionContrasts[,groupInteraction := paste0(names(treatmentFactorsList)[i], " vs ", names(treatmentFactorsList)[j])]
  df_interactionContrasts <- df_interactionContrasts %>% dplyr::mutate(groupInteraction = paste0(names(treatmentFactorsList)[i], " vs ", names(treatmentFactorsList)[j]))
}

#' define all interaction contrasts
#'
#' @param treatmentFactorsList
#' @param groupInteractionToKeep
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
defineAllInteractionContrasts <- function(treatmentFactorsList, groupInteractionToKeep = NULL){
  allInteractionsContrasts_df <- data.table(contrast = character(), groupComparison = factor(), groupInteraction = character(),
                                            outsideGroup = character(),contrastName = character(), type = character())
  # combn(names(treatmentFactorsList),2)
  cat(paste("\ntreatment factors names:\n"))
  print(as.character(names(treatmentFactorsList)))
  vecFori <- combn(length(treatmentFactorsList),2)[1,]
  vecForj <- combn(length(treatmentFactorsList),2)[2,]
  for(k in seq_along(vecFori)){
    i <- vecFori[k]
    j <- vecForj[k]
    #print(paste("i:",i))
    #print(paste("j:",j))
    dataTableToCreate <- defineInteractionConstrastForPairsOfFactors(treatmentFactorsList, i, j)
    allInteractionsContrasts_df <- rbind(allInteractionsContrasts_df, dataTableToCreate)
  }
  if(!missing(groupInteractionToKeep)){
    allInteractionsContrasts_df <- subset(allInteractionsContrasts_df, (groupInteraction %in% groupInteractionToKeep))
  }
  return(allInteractionsContrasts_df)
}


################## function for getContrastMatrix ExpDesign method  ########################


#' compute group binary vector
#'
#' Compute a binary vector (a single vector of 0 and 1) returning the matched group string(s) from a grepl match on the design model matrix colnames.
#'
#' @param biologicalGroups: factor giving group membership (treatment condition associated to one sample)
#' @param colnamesMatrixDesign : vector giving the column names of the model design matrix
#' @param interactionPresent: logical. If TRUE interaction is include in the design model matrix
#'
#' @return a binary vector (a single vector of 0 and 1) returning the matched group string(s) from a grepl match on the design model matrix colnames
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
computeGroupVector <- function(treatmentGroups, colnamesMatrixDesign, interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction) {
  vectorLength <- length(colnamesMatrixDesign)
  groupVector <- rep(0, vectorLength)
  if(interactionPresent){ 
    simples <- unique(unlist(strsplit(treatmentGroups, "_")))
    order2interaction <- combn(simples,2, FUN=paste, collapse=':')
    toMatchList <- c(simples, order2interaction)
    if(isThreeOrderInteraction){
      order3interaction <- combn(simples,3, FUN=paste, collapse=':')
      toMatchList <- c(toMatchList, order3interaction)
    } 
    toMatchList <- paste("^", toMatchList, "$", sep = "")
    #grepl return TRUE for each matched pattern
    pos <- grepl(paste(toMatchList, collapse = "|"), x= colnamesMatrixDesign)
    groupVector <- as.numeric(pos)
    # for intercept
    groupVector[1] <- 1
  } else {
    simples <- unlist(strsplit(treatmentGroups, "_"))
    toMatchList <- simples
    toMatchList <- paste("^", toMatchList, "$", sep = "")
    #toMatch <- paste("^", treatmentGroups, "$", sep = "")
    pos <- grepl(toMatchList, x= colnamesMatrixDesign)
    groupVector <- as.numeric(pos)
    # for intercept
    groupVector[1] <- 1
  }
  return(groupVector)
}
#' Assign binary vector to groups
#'
#' @param treatmentFactorsList: list of treatment factors levels
#' @param interaction_in_model: logical. If TRUE interaction is include in the design model matrix
#' @param modelMatrix: numeric matrix giving the design matrix of the GLM.
#' @param treatmentCondenv: the environment to use
#'
#' @return
#' @export
#'
#' @examples
#' @author Christine Paysant-Le Roux
assignVectorToGroups <- function(treatmentFactorsList = treatmentFactorsList, modelMatrix = modelMatrix,  interactionPresent = interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction, treatmentCondenv = treatmentCondenv){
  # treatment conditions (group) compatible with colnames of the design model matrix
  treatmentGroups <- do.call(paste, c(expand.grid(treatmentFactorsList), sep = "_"))
  # assign binary vector to each group
  modelMatrixColnames <- colnames(modelMatrix)
  # treatmentCondenv <- new.env()
  binaryVectorsList <- lapply(treatmentGroups, function(x) computeGroupVector(x, modelMatrixColnames, interactionPresent, isThreeOrderInteraction = isThreeOrderInteraction))
  names(binaryVectorsList) <- treatmentGroups
  groupDF <- as.data.frame(binaryVectorsList, row.names = modelMatrixColnames)
  mapply(function(x, value) assign(x, value, pos = treatmentCondenv), treatmentGroups, binaryVectorsList)
  # http://r.789695.n4.nabble.com/Using-assign-with-mapply-td4681789.html
  # for (i in 1:n) assign(levels[i], indicatorForModelWithInteraction(levels[i], colnamesGLMdesign = colnamesGLMdesign), pos = levelsenv)
}

#' return contrast coefficients
#'
#' @param contrast
#' @param colnamesGLMdesign
#' @param treatmentCondenv: the environment to use
#' 
#' @return
#' @export
#' @examples
#' @author Christine Paysant-Le Roux
returnContrastCoefficients <- function(contrast, colnamesGLMdesign, treatmentCondenv){
  expression <- NULL
  if (!is.null(contrast)) {
    expression <- as.character(contrast)
    contrastVector <-rep(0,length(colnamesGLMdesign))
    ej <- parse(text = expression[1])
    # eval evaluates the expression argument in the environment specified by envir and returns the computed value
    contrastVector <- eval(ej, envir = treatmentCondenv)
  }
  return(contrastVector)
}

# convenience function to make more descriptive names for factor variables
# function from DESeq2 package
renameModelMatrixColumns <- function(sampleData, modelFormula) {
  data <- as.data.frame(sampleData)
  allVars <- all.vars(modelFormula)
  predictorVars <- allVars[(allVars %in% labels(terms(modelFormula)))]
  colNamesFrom <- make.names(do.call(c,lapply(predictorVars, function(v) paste0(v,levels(data[[v]])[-1]))))
  colNamesTo <- make.names(do.call(c,lapply(predictorVars, function(v) paste0(v,"_",levels(data[[v]])[-1],"_vs_",levels(data[[v]])[1]))))
  data.frame(from=colNamesFrom,to=colNamesTo,stringsAsFactors=FALSE)
}


####################################### CO-EXPRESSION ##############################

#' return gene list
#'
#' @param matrix
#' @param colnames
#' @param mergeType
#' @return list of genes
#' @export
#' @examples
getDEGlist_for_coseqAnalysis <- function(matrix, colnames = colnames(matrix)[-1], mergeType="union"){
  
  if (length(colnames) == 0 ){ return(NULL) }
  
  matrix_sum <- matrix %>% dplyr::mutate(sum = dplyr::select(., all_of(colnames)) %>% rowSums(.))
  
  DEG_list <- switch(mergeType,
                     
         "union"={        dplyr::filter(matrix_sum, sum != 0)[1] },
         "intersection"={ dplyr::filter(matrix_sum, sum == length(colnames))[1] }
  )
  
  if (length(DEG_list$DEG) == 0 ){ return(NULL) }
  
  return(DEG_list$DEG)
}


#' @title run Coseq for co-expression analysis
#' @param counts matrix
#' @param K
#' @param iter
#' @param model 
#' @param transformation
#' @param normFactors 
#' @return coseqResults
#' @export 
runCoseq <- function(counts, K=2:20, iter = 5, model="Normal", transformation="arcsin",  normFactors="TMM"){
            
  
  
            coseq.res <- coseq::coseq(counts, K=K, iter=iter, model=model, transformation=transformation,
                                      parallel=TRUE, meanFilterCutoff=50, normFactors=normFactors, seed=12345)
            
            # Results.1 <- list()
            # Results.1_min_icl <- list()
            # 
            # for (a in 1:iter){
            #   Results.1[[a]] <- coseq(counts, K=K, model=model, transformation=transformation, 
            #                           parallel=parallel, meanFilterCutoff=meanFilterCutoff, normFactors=normFactors)
            #  
            #   Results.1_min_icl[[a]] <- min(metadata(Results.1[[a]])$ICL,na.rm=TRUE)
            # }
            # 
            # coseq.res <- Results.1[[which.min(Results.1_min_icl)]]
            
            return(coseq.res)
          }



#' @title profile boxplot per cluster
#' @param coseq.res coseq object
#' @param selectedCluster cluster num
#' @param conds condition matrix
#' @return 
#' @export 
coseq.y_profile.one.plot <- function(coseq.res, selectedCluster, conds){
  
  nb_cluster <- coseq.res@metadata$nbCluster[min(coseq.res@metadata$ICL) == coseq.res@metadata$ICL]
  
  y_profiles <- list()
  for (i in 1:nb_cluster){
    y_profiles[[i]] <- coseq.res@y_profiles[coseq.res@allResults[[nb_cluster-1]][,i] != 0,] %>% 
      data.frame() %>% reshape2::melt() %>%  dplyr::rename(samples = variable) %>% 
      full_join(conds , by = "samples") %>% dplyr::mutate(cluster = i)
  }
  y_profiles.gg <-  y_profiles %>% purrr::reduce(rbind)
  y_profiles.gg$groups <- factor(y_profiles.gg$groups, levels = unique(conds$groups))
  y_profiles.gg$samples <- factor(y_profiles.gg$samples, levels = unique(conds$samples))
  
  p <- ggplot(data = dplyr::filter(y_profiles.gg, cluster == selectedCluster)) +
    geom_boxplot(aes(x=samples, y=value, fill = groups), outlier.size = 0.3) + facet_wrap(~cluster) +
    theme(axis.text.x=element_blank())
    #theme(axis.text.x=element_text(angle=90, hjust=1))

  print(p)
}



################################### ANNOTATION #############################

#' @title EnrichmentHyperG
#' @param alpha 
#' @param annotation gene annotation
#' @param geneList gene list 
#' @return list
#' @export EnrichmentHyperG
#'
EnrichmentHyperG <- function(annotation, geneList, alpha = 0.01){
#enrichment_analysis <- function(Reference,Gene_List,Alpha){
    
  # ## success in the urn /  For each annotation term, number of annotated genes in the Reference file
  # m=table(Reference[,2])
  # ## failures in the urn / For each annotation term, number of not annotated genes in the Reference file
  # n=length(unique(Reference[,1]))-m
  # 
  # trial<-merge(Gene_List,Reference)
  # ## trial effective / number of genes in the gene list
  # k=length(unique(trial[,1]))
  # ## trial success /  For each annotation term, number of annotated genes in the gene list file
  # x=table(factor(trial[,2],levels=rownames(m)))
  
  ## success in the urn /  For each annotation term, number of annotated genes in the Reference file
  Urn_Success <- annotation %>% dplyr::group_by(Term, Name, Domain) %>% count(name = "Urn_Success")
  
  ## size of reference / nbr of genes in Ref file
  Urn_effective <- length(unique(annotation$geneID))
  
  ## 
  #trial<-merge(geneList,annotation)
  trial <- dplyr::filter(annotation , geneID %in% geneList)
  
  ## trial effective / number of genes in the gene list
  Trial_effective <- length(unique(trial$geneID))
  
  ## trial success /  For each annotation term, number of annotated genes in the gene list file
  Trial_Success <- trial %>% dplyr::group_by(Term, Name, Domain) %>% count(name = "Trial_Success")
            
  
  ## Result files
  res=NaN
  # Term=rownames(m)
  # m=as.numeric(m)
  # n=as.numeric(n)
  # x=as.numeric(x)
  # res=data.frame(Term,Urn_Success=m,Urn_Failures=n,Trial_Success=x,Trial_effective=k,
  #                Urn_percentage_Success=signif(100*m/(m+n),3),
  #                Trial_percentage_Success=signif(100*x/k,3), 
  #                Pvalue_over=phyper(x-1,m,n,k,lower.tail=FALSE),
  #                Pvalue_under=phyper(x,m,n,k,lower.tail=TRUE))       
  
  res= dplyr::full_join(Urn_Success, Trial_Success, by = c("Term", "Name", "Domain")) %>% 
       dplyr::mutate(Urn_percentage_Success   = signif(100*Urn_Success/Urn_effective, 3), 
                     Trial_percentage_Success = signif(100*Trial_Success/Trial_effective, 3), 
                     Pvalue_over  = phyper(Trial_Success-1,Urn_Success, (Urn_effective-Urn_Success),Trial_effective,lower.tail=FALSE), 
                     Pvalue_under = phyper(Trial_Success,  Urn_Success, (Urn_effective-Urn_Success),Trial_effective,lower.tail=TRUE))
  
  # res_over_under <-NULL
  # index=which(res$Pvalue_over<Alpha)
  # if(length(index)!=0)
  # {
  #   res_over <- res[index,]
  #   res_over[,10] <- "overrepresented"
  #   colnames(res_over)[10] <- c("Decision")
  #   res_over_under <- res_over
  # }
  # 
  # index=which(res$Pvalue_under<Alpha)
  # if(length(index)!=0)
  # {
  #   res_under <- res[index,]
  #   res_under[,10] <- "underrepresented"
  #   colnames(res_under)[10] <- c("Decision")
  #   res_over_under <- rbind(res_over_under,res_under)
  # }
  
  res_over_under <- NULL
  res_over_under <- res %>% dplyr::mutate(Decision = dplyr::if_else(Pvalue_over <alpha, "overrepresented", 
                                                     dplyr::if_else(Pvalue_under<alpha, "underrepresented", NULL))) %>%
                            dplyr::filter(!is.na(Decision))
  
  Results <- list("All_results"   = res, "Over_Under_Results" = res_over_under, 
                  "Urn_effective" = Urn_effective, "Trial_effective" = Trial_effective)
            
  return(Results)
            
}



#' pvalue.enrichment.plot
#'
#' @param data 
#' @param Over_Under
#' @param index result size index
#' @param pngFile 
#' @return plot
#' @export
#'
#' @examples
pvalue.enrichment.plot <- function(data, Over_Under, index = "Top50" , pngFile=NULL){
  
  data_ord <- switch (Over_Under,  
          "overrepresented"  = { 
                   dplyr::filter(data, Decision == Over_Under) %>%  dplyr::arrange(desc(Pvalue_over)) %>% 
                   dplyr::mutate(Pvalue = Pvalue_over) 
            },
          "underrepresented" = { 
                   dplyr::filter(data, Decision == Over_Under) %>%  dplyr::arrange(desc(Pvalue_under)) %>% 
                   dplyr::mutate(Pvalue = Pvalue_under) 
            }
          )
  data_ord$Term <- factor(data_ord$Term, levels = data_ord$Term)
  
  switch (index,
    "all" = { p <- ggplot(data = data_ord, aes(x=Pvalue, y=Term, size=Trial_Success, color=Domain)) +  
                     geom_point(alpha=0.5) + scale_size(range = c(0.1, 10))},
    "Top50" = { p <- ggplot(data = tail(data_ord, n=50), aes(x=Pvalue, y=Term, size=Trial_Success, color=Domain)) +  
                     geom_point(alpha=0.5) + scale_size(range = c(0.1, 10))}
  )
  
  print(p)
  
  if (! is.null(pngFile)){
    ggsave(filename = pngFile, plot = p)
  }
  
}


