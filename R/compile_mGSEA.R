compile_mGSEA <- function(ssGSEA.list, p=0.05, FDR=0.25, n.tile.sets=10){
  ## create heatmap data frames
  # extract GSEA results for each omics type
  types <- names(ssGSEA.list)
  GSEA.df <- list()
  for (i in 1:length(types)) {
    GSEA.df[[types[i]]] <- ssGSEA.list[[types[i]]]$result
  }
  
  # collapse GSEA results across omics types
  GSEA.df <- data.table::rbindlist(GSEA.df, use.names = TRUE, idcol = "type")
  GSEA.df$minusLogP <- -log(GSEA.df$p_value, base = 10)
  GSEA.df$minusLogFDR <- -log(GSEA.df$FDR_q_value, base = 10)
  
  # reduce plot data down to top results
  sig.GSEA.df <- GSEA.df[GSEA.df$p_value < p & 
                           GSEA.df$FDR_q_value < FDR, ]
  top.sig.GSEA.df <- sig.GSEA.df %>% dplyr::slice_max(abs(NES), n = n.plot)
  top.GSEA.df <- GSEA.df[GSEA.df$Gene_set %in% top.sig.GSEA.df$Gene_set, ]
  
  ## create tile plot
  # set order of drug sets (decreasing by NES)
  top.GSEA.df <- dplyr::arrange(top.GSEA.df, desc(NES))
  
  # set theme
  bg.theme <- ggplot2::theme(legend.background = element_rect(), legend.position="top", 
                             legend.text = element_text(size=14), 
                             legend.key = element_blank(), 
                             legend.title = element_text(size=16),
                             axis.title.x = element_text(size=20), 
                             axis.text.x  = element_text(size=16, colour = "black"),
                             axis.title.y = element_text(size=20), 
                             axis.text.y  = element_text(size=16, colour = "black"),
                             plot.title = element_text(
                               lineheight=.8, face="bold", size=36), 
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(),
                             panel.border = element_rect(fill=NA), 
                             panel.background = element_blank(),
                             axis.line = element_line(colour = "black"), 
                             axis.ticks.x = element_line(colour="black"),
                             axis.ticks.y = element_line(colour="black"))
  
  # generate tile plot
  tile.plot <- ggplot2::ggplot(top.GSEA.df, 
                               ggplot2::aes(x = type, y = Gene_set, color = NES, 
                                            size = -log10(FDR_q_value))) + 
    ggplot2::geom_point() + 
    ggplot2::scale_y_discrete(limits = top.GSEA.df$Gene_set) + 
    viridis::scale_color_viridis() + bg.theme + ggplot2::labs(x = "Omics Type", 
                                                              y = "Drug Mechanism of Action", 
                                                              color = "NES", size = "-log(FDR)")
  
  ## run correlations
  # extract NES, -logFDR values across omics types
  NES.df <- reshape2::dcast(GSEA.df, Gene_set ~ type, 
                            value.var = "NES", fill = NA)
  minusLogFDR.df <- reshape2::dcast(GSEA.df, Gene_set ~ type, 
                                    value.var = "minusLogFDR", fill = NA)
  
  # convert from data frame to numeric matrix
  rownames(NES.df) <- NES.df$Gene_set
  NES.mat <- as.matrix(NES.df[ , 2:ncol(NES.df)])
  
  # create correlation matrix
  corr.mat <- stats::cor(NES.mat)
  
  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  
  ## compile outputs
  outputs <- list(results = GSEA.df,
                  NES.df = NES.df,
                  minusLogFDR.df = minusLogFDR.df,
                  tile.plot = tile.plot,
                  corr = corr.mat,
                  corr.plot = corr.mat.plot)
  return(outputs)
}
