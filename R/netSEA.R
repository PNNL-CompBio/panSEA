netSEA <- function(inputs, outputs, element.name = rep("Gene", length(inputs)), 
                   rank.var = rep("Log2FC", length(inputs)), 
                   n.top.sets = 3, p = 0.05, FDR = 0.25){
  message("Generating network graph...")
  
  #### Step 1. Identify top sets ####
  outputs <- data.table::rbindlist(outputs, use.names = TRUE, 
                                   idcol = "type")
  sig.outputs <- outputs[outputs$p_value < p & outputs$FDR_q_value < FDR, ]
  top.outputs <- sig.outputs %>% dplyr::slice_max(abs(NES), n = n.top.sets)
  
  #### Step 2. Identify leading edge elements and shared sets ####
  pair <- c()
  leads <- c()
  for (i in 1:nrow(top.outputs)) {
    # get leading edge elements from each set
    set.leads <- stringr::str_split(top.outputs$Leading_edge[i], ", ")[[1]]
    leads <- c(leads, set.leads)
    for (j in 1:length(set.leads)) {
      for (k in 1:length(set.leads)) {
        alpha.leads <- sort(c(set.leads[j], set.leads[k]))
        alpha.pair <- paste0(alpha.leads[1], "_&_", alpha.leads[2])
        pair <- c(pairs, alpha.pair)
      }
    }
  }
  pair <- unique(pairs)
  leads <- unique(leads)
  
  # organize shared set info in data frame
  edge.df <- as.data.frame(pairs)
  edge.df$source <- stringr::str_split_i(edge.df$pair, "_&_", 1)
  edge.df$target <- stringr::str_split_i(edge.df$pair, "_&_", 2)
  edge.df$importance <- NA # this value will determine edge width
  
  #### Step 3. Extract data for graph ####
  # compile inputs based on rank.var
  all.inputs <- inputs[[1]][ , c(element.name[1], rank.var[1])]
  colnames(all.inputs) <- c("feature", "rank")
  for (i in 2:length(inputs)) {
    temp.input <- inputs[[i]][ , c(element.name[i], rank.var[i])]
    colnames(temp.inputs) <- c("feature", "rank")
    all.inputs <- rbind(all.inputs, temp.input)
  }
  
  # get data for nodes (leading edge elements)
  lead.inputs <- all.inputs[all.inputs$feature %in% leads, ]
  
  # compile node info
  node.df <- plyr::ddply(lead.inputs, .(feature), summarize,
                         AvgRank = mean(rank))
  
  ## run correlations between nodes for edge thickness if relevant
  if (nrow(lead.wide) >= 3) {
    # create matrix for correlations
    lead.mat <- reshape2::dcast(lead.inputs, type ~ feature)
    rownames(lead.mat) <- lead.mat$type
    lead.mat$type <- NULL
    lead.mat <- as.matrix(lead.mat)
    
    # create correlation matrix
    cor.result <- stats::cor(lead.mat) # check format to change colnames below
    for (i in 1:nrow(edge.df)) {
      edge.df$importance[i] <- cor.result[
        cor.result$var1 == edge.df$source[i] & 
          cor.result$var2 == edge.df$target[i], ]$cor
    }
  } else {
    net.df$importance <- 1
  }
  
  #### Step 3. Generate network graph ####
  # make color palette
  color.pal <- c("blue", "red") # all the avg ranks should be nonzero
  node.df$carac <- NA
  node.df[node.df$AvgRank > 0, ]$carac <- "Positive"
  node.df[node.df$AvgRank < 0, ]$carac <- "Negative"
  levels(node.df$carac) <- c("Positive", "Negative") # blue = positive and so on
  
  # turn data frame into igraph object
  network <- igraph::graph_from_data_frame(
    d = edge.df, directed = FALSE, vertices = node.df[ , c("feature", "carac")])
  
  # assign node size based on degree of connectivity
  igraph::V(network)$size <- igraph::degree(network, mode = "all")*3
  
  # generate plot
  netPlot <- plot(network) + 
    legend("bottomleft", legend = levels(V(network)$carac), 
           col = color.pal, bty = "n", pch = 20, pt.cex = 3, cex = 1.5, 
           text.col = color.pal, horiz = FALSE, inset = c(0.1, 0.1))

  return(netPlot)
}
