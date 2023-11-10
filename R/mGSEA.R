mGSEA <- function(data, gmt=NULL, direction.adjust=NULL, FDR=0.25, 
                  num.permutations=1000, stat.type="Weighted", min.per.set=6, 
                  sep = ", ", 
                  exclusions = c("-666", "NA", "na", "NaN", "NULL"), 
                  descriptions=NULL, heatmap = FALSE){
  #### Step 1. Perform ssGSEA on each omics type ####
  types <- colnames(data)[2:ncol(data)]
  ssGSEA.list <- list()
  for(i in 1:length(types)) {
    message(paste("Running ssGSEA using", types[i], "data"))
    
    ssGSEA.list[[types[i]]] <- lapply(as.list(data[ , i+1]), ssGSEA, 
                                    gmt, direction.adjust, FDR,
                                    num.permutations, stat.type, min.per.set,
                                    sep, exclusions, descriptions)
  }
  #### Step 2. Compile ssGSEA results across omics types ####
  ## create heatmap data frames
  # extract ssGSEA results for each omics type
  ssGSEA.results <- list()
  for (i in 1:length(types)) {
    ssGSEA.results[[types[i]]] <- ssGSEA.list[[types[i]]]$result
  }
  
  # collapse ssGSEA results across omics types
  ssGSEA.df <- data.table::rbindlist(ssGSEA.results, use.names = TRUE, idcol = "type")
  ssGSEA.df$minusLogP <- -log(ssGSEA.df$p_value, base = 10)
  ssGSEA.df$minusLogFDR <- -log(ssGSEA.df$FDR_q_value, base = 10)
  
  # extract NES, -logP, -logFDR values across omics types
  NES.df <- reshape2::dcast(ssGSEA.df, Gene_set ~ type, value.var = "NES", fill = NA)
  minusLogP.df <- reshape2::dcast(ssGSEA.df, Gene_set ~ type, value.var = "minusLogP", fill = NA)
  minusLogFDR.df <- reshape2::dcast(ssGSEA.df, Gene_set ~ type, value.var = "minusLogFDR", fill = NA)
  
  ## create heatmap
  # get min & max color values
  abs.min.NES <- abs(floor(min(ssGSEA.df$NES)))
  abs.max.NES <- abs(ceiling(max(ssGSEA.df$NES)))
  limit.NES <- ifelse(abs.max.NES > abs.min.NES, abs.max.NES, abs.min.NES)
  
  # generate heatmap
  if (heatmap) {
    heatmap <- morpheus::morpheus(NES.df, 
                                  colorScheme = list(
                                    values = list(-limit.NES, 0, limit.NES)))
  }
  
  # ## run correlations
  # # create correlation matrix
  # corr.mat <- stats::cor(NES.df)
  # 
  # # plot correlation matrix
  # corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  # 
  # #### Step 4. Compile drug results across omics types ####
  # ## create heatmap data frames
  # # extract ssGSEA results for each omics type
  # drug.results <- list()
  # for (i in 1:length(types)) {
  #   drug.results[[i]] <- ssGSEA.list[[types[i]]]$corr.result
  # }
  # 
  # # collapse drug results across omics types
  # drug.df <- data.table::rbindlist(drug.results, use.names = TRUE, idcol = "type")
  # drug.df$minusLogP <- -log(drug.df$Pearson.p, base = 10)
  # drug.df$minusLogFDR <- -log(drug.df$Pearson.q, base = 10)
  # 
  # # extract NES, -logP, -logFDR values across omics types
  # est.df <- reshape2::dcast(drug.df, Drug ~ type, value.var = rank.metric, fill = NA)
  # drug.minusLogP.df <- reshape2::dcast(drug.df, Drug ~ type, value.var = "minusLogP", fill = NA)
  # drug.minusLogFDR.df <- reshape2::dcast(drug.df, Drug ~ type, value.var = "minusLogFDR", fill = NA)
  # 
  # ## create heatmap
  # # get min & max color values
  # abs.min.est <- abs(floor(min(drug.df[ , c(rank.metric)])))
  # abs.max.est <- abs(ceiling(max(drug.df[ , c(rank.metric)])))
  # limit.est <- ifelse(abs.max.est > abs.min.est, abs.max.est, abs.min.est)
  # 
  # # generate heatmap
  # drug.heatmap <- morpheus::morpheus(est.df, 
  #                                    colorScheme = list(
  #                                      values = list(-limit.NES, 0, limit.est)))
  # 
  # ## run correlations
  # # create correlation matrix
  # drug.corr.mat <- corrr:cor(est.df)
  # 
  # # plot correlation matrix
  # drug.corr.mat.plot <- ggcorrplot::ggcorrplot(drug.corr.mat)
  
  return(list(all.results = ssGSEA.list,
              results = ssGSEA.df,
              NES.df = NES.df,
              minusLogP.df = minusLogP.df,
              minusLogFDR.df = minusLogFDR.df
              #,
              # heatmap = heatmap,
              # corr = corr.mat,
              # corr.matrix = corr.mat.plot,
              # drug.est.df = est.df,
              # drug.minusLogP.df = minusLogP.df,
              # drug.minusLogFDR.df = minusLogFDR.df,
              # drug.heatmap = drug.heatmap,
              # drug.corr = drug.corr.mat,
              # drug.corr.matrix = drug.corr.mat.plot
              ))
}
