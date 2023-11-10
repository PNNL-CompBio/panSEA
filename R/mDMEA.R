# Multiomic DMEA
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
# Created: 2023-10-16
# Last modified: 2023-11-10

# input: molecular signature(s) of interest, molecular screening data set(s), drug screen data set, drug annotations
# output: DMEA results for each molecular signature, heatmaps, correlation matrices

mDMEA <- function(drug.sensitivity, gmt=NULL, expression, weights, types, value="AUC",
                  sample.names=colnames(drug.sensitivity)[1], gene.names=colnames(weights[[1]])[1],
                  weight.values=colnames(weights[[1]])[2], rank.metric="Pearson.est", FDR=0.25,
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
  
  # check that sample names are in expression data frames
  for(i in 1:length(types)) {
    if (!(sample.names %in% names(expression[[i]]))) {
      stop("sample.names must match across drug.sensitivity and expression data frames")
    }
  }
  
  #### Step 2. Perform DMEA on each omics type ####
  DMEA.list <- list()
  for(i in 1:length(types)) {
    message(paste("Running DMEA using", types[i], "data"))
    
    DMEA.list[[types[i]]] <- DMEA::DMEA(drug.sensitivity, gmt, expression[[i]], weights[[i]],
                                        value, sample.names, gene.names,
                                        weight.values, rank.metric, FDR,
                                        num.permutations, stat.type, drug.info, 
                                        drug, set.type, min.per.set, sep,
                                        exclusions, descriptions, min.per.corr, 
                                        scatter.plots, scatter.plot.type, 
                                        FDR.scatter.plots, xlab, ylab,
                                        position.x, position.y, se)
    
    # merge correlation results with drug annotations if !is.null(drug.info)
    if (!is.null(drug.info)) { 
      DMEA.list[[types[i]]]$corr.result <- 
        merge(DMEA.list[[types[i]]]$corr.result, drug.info, 
              by.x = "Drug", by.y = drug)
    }
  }
  
  #### Step 3. Compile DMEA results across omics types ####
  compiled.DMEA <- compilDMEA(DMEA.list)
  
  return(list(compiled.results = compiled.DMEA, 
              all.results = DMEA.list
              ))
}
