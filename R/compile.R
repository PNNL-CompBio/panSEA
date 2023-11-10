compilEA <- function(DEGs, ssGSEA.results, DMEA.results, 
                     PTM.DEGs, KSEA.results){
  types <- names(DEGs)
  #### Step 1. Compile differential expression results ####
  DEGs <- data.table::rbindlist(DEGs, use.names = TRUE, idcol = "type")
  PTM.DEGs <- data.table::rbindlist(PTM.DEGs, use.names = TRUE, idcol = "type")

  panDEGs <- plyr::ddply(DEGs, .("Gene"), summarize, 
                         avgLog2FC = mean(LogFC, na.rm = TRUE),
                         sdLog2FC = stats::sd(LogFC, na.rm = TRUE),
                         Fisher_p = as.numeric(metap::sumlog(p-value)$p))
  
  #### Step 2. Compile enrichment results ####
  ssGSEA.df <- list()
  DMEA.df <- list()
  for (i in 1:length(types)) {
    ssGSEA.df[[types[i]]] <- ssGSEA.results[[types[i]]]$result
    DMEA.df[[types[i]]] <- DMEA.results[[types[i]]]$result
  }
  
  
  if (length(types) > 1) {
    ## create heatmap data frames
    # extract DMEA results for each omics type
    DMEA.results <- list()
    for (i in 1:length(types)) {
      DMEA.results[[types[i]]] <- DMEA.list[[types[i]]]$result
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
    if (heatmap) {
      heatmap <- morpheus::morpheus(NES.df, 
                                    colorScheme = list(
                                      values = list(-limit.NES, 0, limit.NES)))
    }
    
    ## run correlations
    # create correlation matrix
    corr.mat <- stats::cor(NES.df)
    
    # plot correlation matrix
    corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
    
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
    
    ## run correlations
    # create correlation matrix
    drug.corr.mat <- corrr:cor(est.df)
    
    # plot correlation matrix
    drug.corr.mat.plot <- ggcorrplot::ggcorrplot(drug.corr.mat)
    
    ## compile outputs
    outputs <- list(results = DMEA.df,
                    NES.df = NES.df,
                    minusLogP.df = minusLogP.df,
                    minusLogFDR.df = minusLogFDR.df,
                    heatmap = heatmap,
                    corr = corr.mat,
                    corr.matrix = corr.mat.plot,
                    element.est.df = est.df,
                    element.minusLogP.df = minusLogP.df,
                    element.minusLogFDR.df = minusLogFDR.df,
                    element.heatmap = element.heatmap,
                    element.corr = element.corr.mat,
                    element.corr.matrix = element.corr.mat.plot)
  }
  
  
  return(outputs)
}
