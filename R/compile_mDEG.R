compile_mDEG <- function(DEGs, p = 0.05, FDR.features = 0.05, 
                         n.dot.features = 10) {
  ## create heatmap data frames
  # extract feature names for each omics type
  types <- names(DEGs)
  for (i in 1:length(types)) {
    DEGs[[types[i]]]$feature_name <- colnames(DEGs[[types[i]]])[[1]]
    colnames(DEGs[[types[i]]])[[1]] <- "feature"
  }

  # collapse DEG results across omics types
  DEG.df <- data.table::rbindlist(DEGs, use.names = TRUE, idcol = "type")
  DEG.df$minusLogP <- -log(DEG.df$P.Value, base = 10)
  DEG.df$minusLogFDR <- -log(DEG.df$adj.P.Val, base = 10)
  DEG.df$sig <- FALSE
  DEG.df[DEG.df$P.Value < p & DEG.df$adj.P.Val < FDR.features, ]$sig <- TRUE

  # summarize results for each feature
  mean.DEG.df <- plyr::ddply(DEG.df, .(feature), summarize,
                             mean_Log2FC = mean(Log2FC),
                             sd_Log2FC = sd(Log2FC),
                             Fisher_p = metap::sumlog(na.omit(P.Value))$p,
                             types = paste0(type, collapse = ", "),
                             N_types = length(unique(type)),
                             N_sig = length(sig[sig]),
                             sig_types = paste0(type[sig], collapse = ", "))
  if (length(unique(mean.DEG.df$Fisher_p)) > 1) {
    mean.DEG.df$adj_Fisher_p <- 
      qvalue::qvalue(mean.DEG.df$Fisher_p, pi0=1)$qvalues
  } else {
    mean.DEG.df$adj_Fisher_p <- NA
  }
  
  # order results by Log2FC for dot plot
  mean.DEG.df <- dplyr::arrange(mean.DEG.df, desc(mean_Log2FC))
  
  # reduce plot data down to top results
  sig.DEG.df <- mean.DEG.df[mean.DEG.df$N_sig > 0, ]
  
  if (nrow(sig.DEG.df) >= n.dot.features) {
    top.DEG.df <- 
      sig.DEG.df %>% dplyr::slice_max(abs(mean_Log2FC), n = n.dot.features)
  } else {
    top.DEG.df <- 
      mean.DEG.df %>% dplyr::slice_max(abs(mean_Log2FC), n = n.dot.features)
  }
  dot.df <- DEG.df[DEG.df$feature %in% top.DEG.df$feature, ]

  ## create venn diagram
  # compile significant results for each type in list
  venn.list <- list()
  for (i in 1:length(types)) {
    venn.list[[types[i]]] <- DEG.df[DEG.df$type == types[i] &
                                       DEG.df$sig, ]$feature
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
      x = type, y = feature, color = Log2FC,
      size = -log10(adj.P.Val)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = mean.DEG.df[
      mean.DEG.df$feature %in% top.DEG.df$feature, ]$feature) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      x = "Input Data",
      y = "Feature Set",
      color = "Log2FC", size = "-log(FDR)"
    )

  ## run correlations
  # extract Log2FC, -logFDR values across omics types
  Log2FC.df <- reshape2::dcast(DEG.df, feature ~ type,
    value.var = "Log2FC", fill = NA
  )
  minusLogFDR.df <- reshape2::dcast(DEG.df, feature ~ type,
    value.var = "minusLogFDR", fill = NA
  )

  # convert from data frame to numeric matrix
  rownames(Log2FC.df) <- Log2FC.df$feature
  Log2FC.mat <- as.matrix(Log2FC.df[, 2:ncol(Log2FC.df)])

  # create correlation matrix
  corr.mat <- stats::cor(Log2FC.mat)

  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)

  ## compile outputs
  outputs <- list(
    results = DEG.df,
    mean.results = mean.DEG.df,
    Log2FC.df = Log2FC.df,
    minusLogFDR.df = minusLogFDR.df,
    venn.diagram = venn.plot,
    dot.plot = dot.plot,
    corr = corr.mat,
    corr.matrix = corr.mat.plot
  )
  return(outputs)
}
