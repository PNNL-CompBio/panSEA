netSEA <- function(inputs, outputs,
                   element.names = rep("Gene", length(inputs)),
                   rank.var = rep("Log2FC", length(inputs)),
                   p = 0.05, FDR = 0.25, n.network.sets = 2 * length(inputs),
                   scale = 5) {
  message("Generating network graph...")

  #### Step 1. Identify significantly enriched sets across all types ####
  outputs <- data.table::rbindlist(outputs,
    use.names = TRUE, idcol = "type", fill = TRUE
  )
  sig.outputs <- outputs[outputs$p_value < p & outputs$FDR_q_value < FDR, ]

  #### Step 2. Extract data for & graph leading edge elements of top sets ####
  if (nrow(sig.outputs) == 0) {
    warning("No enrichments passed significance thresholds for network graph")
    netPlot <- NA
    intPlot <- NA
  } else {
    # identify top significantly enriched sets
    if (nrow(sig.outputs) > n.network.sets) {
      top.outputs <- sig.outputs %>% dplyr::slice_max(abs(NES),
        n = n.network.sets
      )
    } else {
      top.outputs <- sig.outputs
    }

    pairs <- c()
    leads <- c()
    for (i in 1:nrow(top.outputs)) {
      # get leading edge elements from each set
      set.leads <- stringr::str_split(top.outputs$Leading_edge[i], ", ")[[1]]
      leads <- unique(c(leads, set.leads))
      for (j in 1:length(set.leads)) {
        for (k in 1:length(set.leads)) {
          # generate list of unique alphabetically-ordered pairs
          alpha.leads <- sort(c(set.leads[j], set.leads[k]))
          alpha.pair <- paste0(alpha.leads[1], "_&_", alpha.leads[2])
          pairs <- unique(c(pairs, alpha.pair))
        }
      }
    }

    # organize shared set info in data frame
    edge.df <- as.data.frame(pairs)
    edge.df$source <- stringr::str_split_i(edge.df$pair, "_&_", 1)
    edge.df$target <- stringr::str_split_i(edge.df$pair, "_&_", 2)
    edge.df <- edge.df[edge.df$source != edge.df$target, ] # remove self pairs
    edge.df$importance <- NA # this value will determine edge width

    # add Element, Rank columns based on element.names, rank.var
    for (i in 1:length(types)) {
      inputs[[i]][, c("Element", "Rank")] <-
        inputs[[i]][, c(element.names[i], rank.var[i])]

      # keep track of which element.names, rank.var were extracted
      inputs[[i]]$Element.name <- element.names[i]
      inputs[[i]]$Rank.var <- rank.var[i]
    }

    # compile inputs based on rank.var
    names(inputs) <- types
    inputs <- data.table::rbindlist(inputs,
      use.names = TRUE, idcol = "type", fill = TRUE
    )

    # get data for nodes (leading edge elements)
    lead.inputs <- inputs[inputs$Element %in% leads, ]

    # compile node info
    node.df <- na.omit(plyr::ddply(lead.inputs, .(Element), summarize,
      AvgRank = mean(Rank)
    ))
    
    # revise lead inputs in case any were removed with na.omit
    # to make sure all vertex names in edge list are in vertex data frame
    lead.inputs <- lead.inputs[lead.inputs$Element %in% node.df$Element, ]

    ## run correlations between nodes for edge thickness if relevant
    if (length(unique(lead.inputs$type)) >= 3) {
      # create matrix for correlations
      lead.mat <- reshape2::dcast(lead.inputs, type ~ Element,
        value.var = "Rank"
      )
      rownames(lead.mat) <- lead.mat$type
      lead.mat$type <- NULL
      lead.mat <- as.matrix(lead.mat)

      # create correlation matrix
      cor.result <- as.data.frame(stats::cor(lead.mat))
      cor.result$source <- rownames(cor.result)
      cor.result.long <- dplyr::distinct(
        reshape2::melt(cor.result,
          id.vars = "source",
          variable.name = "target", value.name = "importance"
        )
      )
      edge.df <- merge(edge.df[, 1:3], cor.result.long, all.x = TRUE)
    } else {
      edge.df$importance <- 1
    }

    # make color palette where blue = positive, red = negative mean rank
    color.pal <- c("blue", "red") # all the avg ranks should be nonzero
    node.df$carac <- NA
    if (nrow(node.df[node.df$AvgRank > 0, ]) > 0) {
      node.df[node.df$AvgRank > 0, ]$carac <- "Positive"
    }
    if (nrow(node.df[node.df$AvgRank < 0, ]) > 0) {
      node.df[node.df$AvgRank < 0, ]$carac <- "Negative"
    }
    levels(node.df$carac) <- c("Positive", "Negative")

    # turn data frame into igraph object
    network <- igraph::graph_from_data_frame(
      d = edge.df[, c("source", "target", "importance")], directed = FALSE,
      vertices = node.df[, c("Element", "carac")]
    )

    # assign node size, color based on degree of connectivity, mean rank
    igraph::V(network)$size <- abs(scale)*abs(node.df$AvgRank)
    igraph::V(network)$color <- ifelse(node.df$carac == "Positive",
      "blue", "red"
    )

    # generate static plot
    igraph::plot.igraph(network, edge.width = abs(edge.df$importance))
    vertex.value.min <- abs(scale)*floor(min(node.df$AvgRank))
    vertex.value.max <- abs(scale)*ceiling(max(node.df$AvgRank))
    vertex.value.abs.max <- max(abs(vertex.value.min), abs(vertex.value.max))
    vertex.value.neg.range <- seq(-vertex.value.abs.max, 0, length.out = 3)
    vertex.value.pos.range <- seq(0, vertex.value.abs.max, length.out = 3)
    vertex.value.range <- unique(c(
      vertex.value.neg.range,
      vertex.value.pos.range
    ))
    vertex.color.range <- c(rep("red", 2), "white", rep("blue", 2))
    plot.rank.name <- ifelse(length(unique(element.names)) == 1,
      rank.var[1], "Rank"
    )
    legend.title <- paste("Mean", plot.rank.name)
    legend("topleft",
      legend = vertex.value.range,
      pt.cex = abs(vertex.value.range) / 5,
      pch = rep(21, length(vertex.value.range)),
      pt.bg = vertex.color.range, title = legend.title
    )
    netPlot <- recordPlot()

    # generate interactive plot
    vis.legend <- data.frame(
      label = c(legend.title, "> 0", "< 0"),
      color.background = c("white", "blue", "red"),
      color.border = c("white", "black", "black"),
      shape = c("dot", "dot", "dot")
    )
    intPlot <- visNetwork::visIgraph(network) %>%
      visNetwork::visEdges(color = list(color = "grey", highlight = "black")) %>%
      visNetwork::visLegend(addNodes = vis.legend, useGroups = FALSE)
  }
  return(list(static = netPlot, interactive = intPlot))
}
