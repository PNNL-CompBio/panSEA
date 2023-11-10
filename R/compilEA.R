compilDMEA <- function(DMEA.results){
  ## create heatmap data frames
  # extract DMEA results for each omics type
  types <- names(DMEA.results)
  DMEA.df <- list()
  for (i in 1:length(types)) {
    DMEA.df[[types[i]]] <- DMEA.results[[types[i]]]$result
  }

  # collapse DMEA results across omics types
  DMEA.df <- data.table::rbindlist(DMEA.df, use.names = TRUE, idcol = "type")
  DMEA.df$minusLogP <- -log(DMEA.df$p_value, base = 10)
  DMEA.df$minusLogFDR <- -log(DMEA.df$FDR_q_value, base = 10)
  
  # extract NES, -logP, -logFDR values across omics types
  NES.df <- reshape2::dcast(DMEA.df, Drug_set ~ type, 
                            value.var = "NES", fill = NA)
  minusLogP.df <- reshape2::dcast(DMEA.df, Drug_set ~ type, 
                                  value.var = "minusLogP", fill = NA)
  minusLogFDR.df <- reshape2::dcast(DMEA.df, Drug_set ~ type, 
                                    value.var = "minusLogFDR", fill = NA)
  
  ## create heatmap
  # get min & max color values
  abs.min.NES <- abs(floor(min(DMEA.df$NES)))
  abs.max.NES <- abs(ceiling(max(DMEA.df$NES)))
  limit.NES <- ifelse(abs.max.NES > abs.min.NES, abs.max.NES, abs.min.NES)
  
  # generate heatmap
  heatmap <- morpheus::morpheus(NES.df, colorScheme = list(
    values = list(-limit.NES, 0, limit.NES)))
  
  ## run correlations
  # create correlation matrix
  corr.mat <- stats::cor(NES.df)
  
  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  
  ## compile outputs
  outputs <- list(results = DMEA.df,
                  NES.df = NES.df,
                  minusLogP.df = minusLogP.df,
                  minusLogFDR.df = minusLogFDR.df,
                  heatmap = hDMEAtmap,
                  corr = corr.mat,
                  corr.matrix = corr.mat.plot)
  return(outputs)
}
