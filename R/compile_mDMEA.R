compile_mDMEA <- function(mDMEA.results, p = 0.05, FDR = 0.25,
                          n.dot.sets = 10) {
  ## compile data
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
  DMEA.df$sig <- FALSE
  DMEA.df[DMEA.df$p_value < p & DMEA.df$FDR_q_value < FDR, ]$sig <- TRUE

  # summarize results for each Drug_set
  mean.DMEA.df <- plyr::ddply(DMEA.df, .(Drug_set), summarize,
                              mean_NES = mean(NES),
                              sd_NES = sd(NES),
                              Fisher_p = as.numeric(metap::sumlog(na.omit(p_value))$p),
                              types = paste0(type, collapse = ", "),
                              N_types = length(unique(type)),
                              N_sig = length(sig[sig]),
                              sig_types = paste0(type[sig], collapse = ", "))
  if (length(unique(mean.DMEA.df$Fisher_p)) > 1) {
    mean.DMEA.df$adj_Fisher_p <- 
      qvalue::qvalue(mean.DMEA.df$Fisher_p, pi0=1)$qvalues
  } else {
    mean.DMEA.df$adj_Fisher_p <- NA
  }
  
  # order results by NES for dot plot
  mean.DMEA.df <- dplyr::arrange(mean.DMEA.df, desc(mean_NES))
  
  # reduce plot data down to top results
  sig.DMEA.df <- mean.DMEA.df[mean.DMEA.df$N_sig > 0, ]
  
  if (nrow(sig.DMEA.df) >= n.dot.sets) {
    top.DMEA.df <- 
      sig.DMEA.df %>% dplyr::slice_max(abs(mean_NES), n = n.dot.sets)
  } else {
    top.DMEA.df <- 
      mean.DMEA.df %>% dplyr::slice_max(abs(mean_NES), n = n.dot.sets)
  }
  dot.df <- DMEA.df[DMEA.df$Drug_set %in% top.DMEA.df$Drug_set, ]

  ## create venn diagram
  # compile significant results for each type in list
  venn.list <- list()
  for (i in 1:length(types)) {
    venn.list[[types[i]]] <- DMEA.df[DMEA.df$type == types[i] &
                                       DMEA.df$sig, ]$Drug_set
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
      x = type, y = Drug_set, color = NES,
      size = -log10(FDR_q_value)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = mean.DMEA.df[
      mean.DMEA.df$Drug_set %in% top.DMEA.df$Drug_set, ]$Drug_set) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      x = "Input Data",
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
    venn.diagram = venn.plot,
    dot.plot = dot.plot,
    corr = corr.mat,
    corr.matrix = corr.mat.plot
  )
  return(outputs)
}
