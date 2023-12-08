compile_mDMEA <- function(mDMEA.results, p = 0.05, FDR = 0.25,
                          n.dot.sets = 10) {
  ## create heatmap data frames
  # extract DMEA results for each omics type
  types <- names(mDMEA.results)
  DMEA.df <- list()
  for (i in 1:length(types)) {
    DMEA.df[[types[i]]] <- mDMEA.results[[types[i]]]$result
  }

  # collapse DMEA results across omics types
  DMEA.df <- data.table::rbindlist(DMEA.df, use.names = TRUE, idcol = "type")
  DMEA.df$minusLogP <- -log(DMEA.df$p_value, base = 10)
  DMEA.df$minusLogFDR <- -log(DMEA.df$FDR_q_value, base = 10)

  # reduce plot data down to top results
  sig.DMEA.df <- DMEA.df[DMEA.df$p_value < p &
    DMEA.df$FDR_q_value < FDR, ]
  
  if (nrow(sig.DMEA.df) > 0) {
    top.sig.DMEA.df <- 
      sig.DMEA.df %>% dplyr::slice_max(abs(NES), n = n.dot.sets)
    top.DMEA.df <- 
      DMEA.df[DMEA.df$Drug_set %in% top.sig.DMEA.df$Drug_set, ]
  } else {
    top.DMEA.df <- 
      DMEA.df %>% dplyr::slice_max(abs(NES), n = n.dot.sets)
  }

  ## create dot plot
  # set order of drug sets (decreasing by mean NES)
  mean.DMEA.df <- plyr::ddply(top.DMEA.df, .(Drug_set), summarize,
                              mean_NES = mean(NES),
                              Fisher_p = as.numeric(metap::sumlog(p_value)$p),
                              types = paste0(type, collapse = ", "),
                              N_types = length(unique(type)))
  mean.DMEA.df$adj_Fisher_p <- qvalue::qvalue(mean.DMEA.df$Fisher_p, pi0=1)$qvalues
  mean.DMEA.df <- dplyr::arrange(
    mean.DMEA.df[mean.DMEA.df$Drug_set %in% top.DMEA.df$Drug_set, ], 
    desc(mean_NES))

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
    top.DMEA.df,
    ggplot2::aes(
      x = type, y = Drug_set, color = NES,
      size = -log10(FDR_q_value)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = mean.DMEA.df$Drug_set) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      x = "Omics Type",
      y = "Drug Mechanism of Action",
      color = "NES", size = "-log(FDR)"
    )

  ## run correlations
  # extract NES, -logFDR values across omics types
  NES.df <- reshape2::dcast(DMEA.df, Drug_set ~ type,
    value.var = "NES", fill = NA
  )
  minusLogFDR.df <- reshape2::dcast(DMEA.df, Drug_set ~ type,
    value.var = "minusLogFDR", fill = NA
  )

  # convert from data frame to numeric matrix
  rownames(NES.df) <- NES.df$Drug_set
  NES.mat <- as.matrix(NES.df[, 2:ncol(NES.df)])

  # create correlation matrix
  corr.mat <- stats::cor(NES.mat)

  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)

  ## compile outputs
  outputs <- list(
    results = DMEA.df,
    mean.results = mean.DMEA.df,
    NES.df = NES.df,
    minusLogFDR.df = minusLogFDR.df,
    dot.plot = dot.plot,
    corr = corr.mat,
    corr.matrix = corr.mat.plot
  )
  return(outputs)
}
