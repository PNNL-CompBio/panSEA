panSEA <- function(data, types, feature.names = rep("Gene", length(types)), 
                   GSEA.rank.var = rep("Log2FC", length(types)), 
                   DMEA.rank.var = rep("Pearson.est", length(types)),
                   group.names=c("Diseased", "Healthy"), 
                   group.samples=c(1:(0.5*ncol(data[[1]])), 
                                   (0.5*ncol(data[[1]])+1):ncol(data[[1]])), 
                   gmt.features=rep("msigdb_Homo sapiens_C2_CP:KEGG", 
                                    length(types)), 
                   gmt.drugs="PRISM", p=0.05, FDR=0.25, num.permutations=1000, 
                   stat.type="Weighted", min.per.set=6, scatter.plots=TRUE,
                   scatter.plot.type="pearson", drug.sensitivity="PRISM", 
                   expression=as.list(rep("adherent CCLE", length(types))),
                   n.network.sets = length(types), n.dot.sets = 10) {
  #### Step 1. Check if formats are correct ####
  # make sure that there are as many types as other inputs
  if (length(types) != length(gmt.features) | 
             length(types) != length(feature.names) | 
             length(types) != length(GSEA.rank.var) | 
             length(types) != length(DMEA.rank.var) | 
             length(types) != length(expression)) {
    stop(paste("Lengths of types, feature.names, GSEA.rank.var,", 
               "DMEA.rank.var, gmt.features, and expression must all match"))
  }
  
  if (length(group.names) > 2 | length(group.names) < 1) {
    stop("Only 1 or 2 group.names are allowed")
  } else if (length(group.names) == 2) {
    #### Step 2. Differential expression analysis ####
    deg <- mDEG(data, types, group.names, group.samples)
    
    #### Step 3. Enrichment analyses ####
    ### ssGSEA & network graph
    if (!is.null(gmt.features)) {
      ssGSEA.results <- mGSEA(deg$DEGs, gmt.features, types, feature.names, 
                              p = p, FDR = FDR, 
                              num.permutations = num.permutations, 
                              stat.type = stat.type, min.per.set = min.per.set, 
                              n.dot.sets = n.dot.sets)
      
      # compile inputs & results for network graph
      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- ssGSEA.results[[types[i]]]$result
      }
      
      ssGSEA.network <- netSEA(deg$DEGs, outputs, feature.names, 
                              GSEA.rank.var, p, FDR, n.network.sets)
    }
    
    ### DMEA & network graph
    if (!is.null(drug.sensitivity) & !is.null(expression) 
        & !is.null(gmt.drugs)) {
      DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, expression, deg, 
                            types, rank.metric = DMEA.rank.var, 
                            weight.values = GSEA.rank.var, p=p, FDR = FDR, 
                            num.permutations = num.permutations, 
                            stat.type = stat.type, min.per.set = min.per.set, 
                            scatter.plots = scatter.plots, 
                            scatter.plot.type = scatter.plot.type, 
                            n.dot.sets = n.dot.sets)
      
      # compile inputs & outputs for network graph
      inputs <- list()
      for (i in 1:length(types)) {
        inputs[[types[i]]] <- DMEA.results[[types[i]]]$corr.result
      }
      
      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- DMEA.results[[types[i]]]$result
      }
      
      DMEA.network <- netSEA(inputs, outputs, rep("Drug", length(inputs)), 
                            DMEA.rank.var, p, FDR, n.network.sets)
    } 
    DEGs <- deg$DEGs
    Log2Transformed <- deg$Log2Transformed
  } else if (length(group.names) == 1) {
    DEGs <- NULL
    Log2Transformed <- NULL
    ssGSEA.results <- list()
    ssGSEA.network <- list()
    DMEA.results <- list()
    DMEA.network <- list()
    for (j in 1:ncol(data)) {
      #### Step 2. Extract omics data for each sample ####
      temp.data <- list()
      for (i in 1:length(types)) {
        temp.df <- data[[i]][ , c(j)]
        temp.df$feature <- rownames(temp.df)
        temp.data[[types[i]]] <- temp.df[ , c(2, 1)]
      }
      
      #### Step 3. Enrichment analyses ####
      ### ssGSEA & network graph
      if (!is.null(gmt.features)) {
        ssGSEA.results[[colnames(data[[1]])[j]]] <- 
          mGSEA(temp.data, gmt.features, types, p = p, FDR = FDR, 
                num.permutations = num.permutations, 
                stat.type = stat.type, min.per.set = min.per.set, 
                n.dot.sets = n.dot.sets)
        
        # compile inputs & outputs for network graph
        outputs <- list()
        for (i in 1:length(types)) {
          outputs[[types[i]]] <- ssGSEA.results[[types[i]]]$result
        }
        
        ssGSEA.network[[colnames(data[[1]])[j]]] <- 
          netSEA(temp.data, outputs, feature.names, 
                GSEA.rank.var, p, FDR, n.network.sets)
      }
      
      ### DMEA & network graph
      if (!is.null(drug.sensitivity) & !is.null(expression) 
          & !is.null(gmt.drugs)) {
        DMEA.results[[colnames(data[[1]])[j]]] <- 
          mDMEA(drug.sensitivity, gmt.drugs, expression, temp.data, types, 
                rank.metric = DMEA.rank.var, weight.values = GSEA.rank.var, 
                p=p, FDR = FDR, num.permutations = num.permutations, 
                stat.type = stat.type, min.per.set = min.per.set, 
                scatter.plots = scatter.plots, 
                scatter.plot.type = scatter.plot.type, 
                n.dot.sets = n.dot.sets)
        
        # compile inputs & outputs for network graph
        inputs <- list()
        for (i in 1:length(types)) {
          inputs[[types[i]]] <- DMEA.results[[types[i]]]$corr.result
        }
        
        outputs <- list()
        for (i in 1:length(types)) {
          outputs[[types[i]]] <- DMEA.results[[types[i]]]$result
        }
        
        DMEA.network[[colnames(data[[1]])[j]]] <- 
          netSEA(inputs, outputs, rep("Drug", length(inputs)), 
                DMEA.rank.var, p, FDR, n.network.sets)
      }
    }
  }
  
  return(list(DEGs = DEGs, Log2Transformed = Log2Transformed, 
              mGSEA.results = ssGSEA.results, mDMEA.results = DMEA.results, 
              mGSEA.network = ssGSEA.network, mDMEA.network = DMEA.network))
}
