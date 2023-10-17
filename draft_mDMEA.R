# Multiomic DMEA
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
# Created: 2023-10-16
# Last modified: 2023-10-16

# input: molecular signature(s) of interest, molecular screening data set(s), drug screen data set, drug annotations
# output: DMEA results for each molecular signature & overall

mDMEA <- function(drug.sensitivity, gmt=NULL, expression, weights, types, value="AUC",
                  sample.names=colnames(expression)[1], gene.names=colnames(weights)[1],
                  weight.values=colnames(weights)[2], rank.metric="Pearson.est", FDR=0.25,
                  num.permutations=1000, stat.type="Weighted", drug.info=NULL, drug="Drug",
                  set.type="moa", min.per.set=6, sep="[|]",
                  exclusions=c("-666", "NA", "na", "NaN", "NULL"), descriptions=NULL,
                  min.per.corr=3, scatter.plots=TRUE, scatter.plot.type="pearson",
                  FDR.scatter.plots=0.05, xlab="Weighted Voting Score", ylab=value,
                  position.x="min", position.y="min", se=TRUE){
  #### Step 1. Check if formats are correct ####
  # check that there are as many types as expression and weights data frames
  if (length(types) != length(expression) | length(types) != length(weights) | 
      length(weights) != length(expression)) {
    stop("Length of types vector must match those of expression and weights lists")
  }
  
  # check that sample names are in drug.sensitivity data frame
  if (!(sample.names %in% names(drug.sensitivity))) {
    stop("sample.names must match across drug.sensitivity and expression data frames")
  }
  
  # check that sample names are in expression and weights data frames
  for(i in 1:length(types)) {
    if (!(sample.names %in% names(expression[[i]]))) {
      stop("sample.names must match across drug.sensitivity and expression data frames")
    } else if (!(sample.names %in% names(weights[[i]]))) {
      stop("sample.names must match across drug.sensitivity and expression data frames")
    }
  }
  
  #### Step 2. Perform DMEA on each omics type ####
  DMEA.list <- list()
  for(i in 1:length(types)) {
    message(paste("Running DMEA using", types[i], "data"))
    
    DMEA.list[[types[i]]] <- DMEA::DMEA(drug.sensitivity, gmt, expression, weights,
                                        value, sample.names, gene.names,
                                        weight.values, rank.metric, FDR,
                                        num.permutations, stat.type, drug.info, 
                                        drug, set.type, min.per.set, sep,
                                        exclusions, descriptions, min.per.corr, 
                                        scatter.plots, scatter.plot.type, 
                                        FDR.scatter.plots, xlab, ylab,
                                        position.x, position.y, se)
  }
  
  
  #### Step 3. Compile DMEA results across omics types ####
  ## create heatmap data frames
  # extract DMEA results for each omics type
  DMEA.results <- list()
  for (i in 1:length(types)) {
    DMEA.results[[i]] <- DMEA.list[[types[i]]]$result
  }
  
  # collapse DMEA results across omics types
  DMEA.df <- data.table::rbindlist(DMEA.results, use.names = TRUE, idcol = "type")
  DMEA.df$minusLogP <- -log(DMEA.df$p_value, base = 10)
  DMEA.df$minusLogFDR <- -log(DMEA.df$FDR_q_value, base = 10)
  
  # extract NES, -logP, -logFDR values across omics types
  NES.df <- reshape2::dcast(DMEA.df, Drug_set ~ type, value.var = "NES", fill = NA)
  minusLogP.df <- reshape2::dcast(DMEA.df, Drug_set ~ type, value.var = "minusLogP", fill = NA)
  minusLogFDR.df <- reshape2::dcast(DMEA.df, Drug_set ~ type, value.var = "minusLogFDR", fill = NA)
  
  ## create heatmap
  # get min & max color values
  abs.min.NES <- abs(floor(min(DMEA.df$NES)))
  abs.max.NES <- abs(ceiling(max(DMEA.df$NES)))
  limit.NES <- ifelse(abs.max.NES > abs.min.NES, abs.max.NES, abs.min.NES)
  
  # generate heatmap
  heatmap <- morpheus::morpheus(NES.df, 
                                colorScheme = list(
                                  values = list(-limit.NES, 0, limit.NES)))
  
  ## run PCA
  # format input data frame
  pca.df <- NES.df
  rownames(pca.df) <- pca.df$Drug_set
  pca.df$Drug_set <- NULL
  
  # calculate components
  pca <- prcomp(pca.df, scale = TRUE)
  
  # store plot of variance explained by each component
  pca.var <- factoextra::fviz_eig(pca)
  
  # store plot of individuals
  pca.plot <- factoextra::fviz_pca_ind(pca, repel = TRUE)
  
  # create UMAP plot
  
  
  # run correlations
  
  #### Step 4. Compile drug results across omics types ####
  ## create heatmap data frames
  # extract DMEA results for each omics type
  drug.results <- list()
  for (i in 1:length(types)) {
    drug.results[[i]] <- DMEA.list[[types[i]]]$corr.result
  }
  
  # collapse drug results across omics types
  drug.df <- data.table::rbindlist(drug.results, use.names = TRUE, idcol = "type")
  drug.df$minusLogP <- -log(drug.df$Pearson.p, base = 10)
  drug.df$minusLogFDR <- -log(drug.df$Pearson.q, base = 10)
  
  # extract NES, -logP, -logFDR values across omics types
  est.df <- reshape2::dcast(drug.df, Drug ~ type, value.var = rank.metric, fill = NA)
  drug.minusLogP.df <- reshape2::dcast(drug.df, Drug ~ type, value.var = "minusLogP", fill = NA)
  drug.minusLogFDR.df <- reshape2::dcast(drug.df, Drug ~ type, value.var = "minusLogFDR", fill = NA)
  
  ## create heatmap
  # get min & max color values
  abs.min.est <- abs(floor(min(drug.df[ , c(rank.metric)])))
  abs.max.est <- abs(ceiling(max(drug.df[ , c(rank.metric)])))
  limit.est <- ifelse(abs.max.est > abs.min.est, abs.max.est, abs.min.est)
  
  # generate heatmap
  drug.heatmap <- morpheus::morpheus(est.df, 
                                colorScheme = list(
                                  values = list(-limit.NES, 0, limit.est)))
  
  ## run PCA
  # calculate components
  
  
  # create UMAP plot
  
  # run correlations
  
  
  return(list(all.results = DMEA.list,
              results = DMEA.df,
              NES.df = NES.df,
              minusLogP.df = minusLogP.df,
              minusLogFDR.df = minusLogFDR.df,
              heatmap = heatmap,
              pca = pca.plot,
              umap = umap.plot,
              corr = DMEA.corr,
              drug.est.df = est.df,
              drug.minusLogP.df = minusLogP.df,
              drug.minusLogFDR.df = minusLogFDR.df,
              drug.heatmap = drug.heatmap,
              drug.pca = drug.pca.plot,
              drug.umap = drug.umap.plot,
              drug.corr = drug.corr))
}