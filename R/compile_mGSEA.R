compile_mGSEA <- function(ssGSEA.list, p = 0.05, FDR = 0.25, n.dot.sets = 10) {
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
  GSEA.df$sig <- FALSE
  GSEA.df[GSEA.df$p_value < p & GSEA.df$FDR_q_value < FDR, ]$sig <- TRUE
  
  # summarize results for each Feature_set
  mean.GSEA.df <- plyr::ddply(GSEA.df, .(Feature_set), summarize,
                              mean_NES = mean(NES),
                              sd_NES = sd(NES),
                              Fisher_p = as.numeric(metap::sumlog(na.omit(p_value))$p),
                              types = paste0(type, collapse = ", "),
                              N_types = length(unique(type)),
                              N_sig = length(sig[sig]),
                              sig_types = paste0(type[sig], collapse = ", "))
  if (length(unique(mean.GSEA.df$Fisher_p)) > 1) {
    mean.GSEA.df$adj_Fisher_p <- 
      qvalue::qvalue(mean.GSEA.df$Fisher_p, pi0=1)$qvalues
  } else {
    mean.GSEA.df$adj_Fisher_p <- NA
  }
  
  # order results by NES for dot plot
  mean.GSEA.df <- dplyr::arrange(mean.GSEA.df, desc(mean_NES))
  
  
  # reduce plot data down to top results
  sig.GSEA.df <- mean.GSEA.df[mean.GSEA.df$N_sig > 0, ]
  
  if (nrow(sig.GSEA.df) >= n.dot.sets) {
    top.GSEA.df <- 
      sig.GSEA.df %>% dplyr::slice_max(abs(mean_NES), n = n.dot.sets)
  } else {
    top.GSEA.df <- 
      mean.GSEA.df %>% dplyr::slice_max(abs(mean_NES), n = n.dot.sets)
  }
  dot.df <- GSEA.df[GSEA.df$Feature_set %in% top.GSEA.df$Feature_set, ]

  ## create venn diagram
  # compile significant results for each type in list
  venn.list <- list()
  for (i in 1:length(types)) {
    venn.list[[types[i]]] <- GSEA.df[GSEA.df$type == types[i] &
                                       GSEA.df$sig, ]$Feature_set
  }
  
  # generate venn diagram
  venn.plot <- ggvenn::ggvenn(venn.list) # only displays first 4 types
  
  ## create dot plot
  # set theme
  bg.theme <- ggplot2::theme(
    legend.background = element_rect(), legend.position = "top",
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    legend.title = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(
      lineheight = .8, face = "bold", size = 36
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )

  # generate dot plot
  dot.plot <- ggplot2::ggplot(
    dot.df,
    ggplot2::aes(
      x = type, y = Feature_set, color = NES,
      size = -log10(FDR_q_value)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = mean.GSEA.df[
      mean.GSEA.df$Feature_set %in% top.GSEA.df$Feature_set, ]$Feature_set) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      x = "Omics Type",
      y = "Feature Set",
      color = "NES", size = "-log(FDR)"
    )

  ## run correlations
  # extract NES, -logFDR values across omics types
  NES.df <- reshape2::dcast(GSEA.df, Feature_set ~ type,
    value.var = "NES", fill = NA
  )
  minusLogFDR.df <- reshape2::dcast(GSEA.df, Feature_set ~ type,
    value.var = "minusLogFDR", fill = NA
  )

  # convert from data frame to numeric matrix
  rownames(NES.df) <- NES.df$Feature_set
  NES.mat <- as.matrix(NES.df[, 2:ncol(NES.df)])

  # create correlation matrix
  corr.mat <- stats::cor(NES.mat)

  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)

  ## compile outputs
  outputs <- list(
    results = GSEA.df,
    mean.results = mean.GSEA.df,
    NES.df = NES.df,
    minusLogFDR.df = minusLogFDR.df,
    venn.diagram = venn.plot,
    dot.plot = dot.plot,
    corr = corr.mat,
    corr.plot = corr.mat.plot
  )
  return(outputs)
}
