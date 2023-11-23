# Multiomic DMEA
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
# Created: 2023-10-16
# Last modified: 2023-11-10

# input: molecular signature(s) of interest, molecular screening data set(s), drug screen data set, drug annotations
# output: DMEA results for each molecular signature, heatmaps, correlation matrices

mDMEA <- function(drug.sensitivity="PRISM", gmt="PRISM", 
                  expression=as.list(rep("adherent CCLE", length(types))), 
                  weights, types, value="AUC", 
                  sample.names=colnames(drug.sensitivity)[1], 
                  gene.names=rep("Gene", length(types)),
                  weight.values=rep("Log2FC", length(types)), 
                  rank.metric=rep("Pearson.est", length(types)), p=0.05, 
                  FDR=0.25, num.permutations=1000, stat.type="Weighted", 
                  drug.info=NULL, drug="Drug", set.type="moa", min.per.set=6, 
                  sep="[|]", exclusions=c("-666", "NA", "na", "NaN", "NULL"), 
                  descriptions=NULL, min.per.corr=3, scatter.plots=TRUE, 
                  scatter.plot.type="pearson", FDR.scatter.plots=0.05, 
                  xlab="Weighted Voting Score", ylab=value, position.x="min", 
                  position.y="min", se=TRUE, n.dot.sets=10){
  #### Step 1. Check if formats are correct ####
  # get drug.sensitivity, gmt, and expression if PRISM/CCLE
  if (drug.sensitivity == "PRISM") {
    PRISM.AUC <- read.csv(file=paste0("https://raw.github.com/BelindaBGarana/", 
    "DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv")) #481 cell lines
    PRISM.AUC$X <- NULL
  }
  
  if (gmt == "PRISM") {
    gmt <- GSA::GSA.read.gmt(
      file=paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/", 
                  "Inputs/MOA_gmt_file_n6_no_special_chars.gmt"))
  }
  
  if ("adherent CCLE" %in% expression) {
    download.file(
      paste0(
        "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/", 
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"), 
      destfile = 
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
    
    download.file(
      paste0(
        "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/", 
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"), 
      destfile = 
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
    
    RNA.df <- rbind(RNA.first200, RNA.rest)
    
    expression[expression == "adherent CCLE"] <- RNA.df
  }
  
  # check that there are as many types as weights data frames
  if (length(types) != length(weights)) {
    stop("Length of types vector must match that of weights list")
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
    
    DMEA.list[[types[i]]] <- DMEA::DMEA(drug.sensitivity, gmt, expression[[i]], 
                                        weights[[i]], value, sample.names, 
                                        gene.names[i], weight.values[i], 
                                        rank.metric[i], FDR, num.permutations, 
                                        stat.type, drug.info, drug, set.type, 
                                        min.per.set, sep, exclusions, 
                                        descriptions, min.per.corr, 
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
  compiled.DMEA <- compile_mDMEA(DMEA.list, p, FDR, n.dot.sets)
  
  return(list(compiled.results = compiled.DMEA, 
              all.results = DMEA.list
              ))
}
