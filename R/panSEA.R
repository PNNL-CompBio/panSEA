panSEA <- function(data.list, types, feature.names = rep("Gene", length(types)),
                   GSEA.rank.var = rep("Log2FC", length(types)),
                   DMEA.rank.var = rep("Pearson.est", length(types)),
                   group.names = c("Diseased", "Healthy"),
                   group.samples = list(
                     2:(0.5 * (ncol(data.list[[1]]) + 1)),
                     (0.5 * (ncol(data.list[[1]]) + 1) + 1):ncol(data.list[[1]])
                   ), group.names2 = NULL, group.samples2 = NULL,
                   gmt.features = as.list(rep(
                     "msigdb_Homo sapiens_C2_CP:KEGG",
                     length(types)
                   )),
                   gmt.drugs = "PRISM", GSEA = TRUE, DMEA = TRUE,
                   DMEA.type = "WV", p = 0.05, FDR = 0.25, 
                   FDR.features = 0.05, num.permutations = 1000, 
                   stat.type = "Weighted", min.per.set = 6, 
                   scatter.plots = TRUE, scatter.plot.type = "pearson", 
                   drug.sensitivity = "PRISM",
                   expression = as.list(rep("adherent CCLE", length(types))),
                   n.network.sets = 2*length(types), n.dot.sets = 10, scale = 5) {
  #### Step 1. Check if formats are correct ####
  # make sure that there are as many types as other inputs
  if (length(types) != length(gmt.features) |
    length(types) != length(feature.names) |
    length(types) != length(GSEA.rank.var) |
    length(types) != length(DMEA.rank.var) |
    length(types) != length(expression)) {
    stop(paste(
      "Lengths of types, feature.names, GSEA.rank.var,",
      "DMEA.rank.var, gmt.features, and expression must all match"
    ))
  }
  
  # make sure grouping lengths are valid
  if (length(group.names) != length(group.samples)) {
    stop("Lengths of group.names must match that of group.samples")
  } else if (length(group.names) > 2 | length(group.names) < 1) {
    stop("Only 1 or 2 group.names are allowed")
  } else if (length(group.names) == 1) {
    DEGs <- NA
    group1.valid <- TRUE
  } else {
    group1.valid <- TRUE
  }
  
  if (length(group.names2) != length(group.samples2)) {
    stop("Lengths of group.names2 must match that of group.samples2")
  } else if (is.null(group.names2) & is.null(group.samples2)) {
    group2.valid <- TRUE
  } else if (length(group.names2) > 2) {
    stop("Only 1 or 2 group.names2 are allowed")
  } else if (length(group.names2) == 1) {
    group2.valid <- TRUE
  } else {
    group2.valid <- TRUE
  }
  
  # check if DMEA.type is valid
  if (DMEA & !(DMEA.type %in% c("WV", "cell_corr", "gene_corr"))) {
    stop(paste("DMEA.type parameter must be either 'WV' for drugs to be ranked",
    "by correlations with weighted voting scores, 'cell_corr' for correlations",
    "with expression values, or 'gene_corr' for correlations with perturbation",
    "values"))
  }
  
  #### Step 2. Differential expression analysis if 2 groups ####
  if (length(group.names) == 2 & group1.valid & group2.valid) {
    DEGs <- panSEA::mDEG(data.list, types, group.names, group.samples,
                         group.names2, group.samples2,
                         feature.names, p, FDR.features, n.dot.sets)
  }

  #### Step 3. Enrichment analyses ####
  if (length(group.names) == 2) {
    ## ssGSEA & network graph
    if (GSEA & !is.null(gmt.features)) {
      ssGSEA.results <- panSEA::mGSEA(DEGs$all.results, gmt.features, types,
        feature.names, GSEA.rank.var, p = p, FDR = FDR, 
        num.permutations = num.permutations, stat.type = stat.type, 
        min.per.set = min.per.set, n.dot.sets = n.dot.sets
      )

      # compile inputs & results for network graph
      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- ssGSEA.results$all.results[[types[i]]]$result
      }

      ssGSEA.network <- panSEA::netSEA(
        DEGs$all.results, outputs, feature.names,
        GSEA.rank.var, p, FDR, n.network.sets, scale
      )
    } else {
      ssGSEA.results <- NA
      ssGSEA.network <- NA
    }

    ## DMEA & network graph
    if (DMEA & !is.null(drug.sensitivity) & !is.null(expression) &
      !is.null(gmt.drugs)) {
      if (DMEA.type == "WV") {
        DMEA.results <- panSEA::mDMEA(drug.sensitivity, gmt.drugs, expression,
                                      DEGs$all.results, types, 
                                      feature.names = feature.names,
                                      rank.metric = DMEA.rank.var, 
                                      weight.values = GSEA.rank.var, 
                                      p = p, FDR = FDR, 
                                      num.permutations = num.permutations,
                                      stat.type = stat.type, 
                                      min.per.set = min.per.set,
                                      scatter.plots = scatter.plots, 
                                      scatter.plot.type = scatter.plot.type,
                                      n.dot.sets = n.dot.sets
        ) 
      } else if (DMEA.type == "cell_corr") {
        DMEA.results <- 
          panSEA::mDMEA_cell_corr(drug.sensitivity, gmt.drugs, expression, 
                                  DEGs$all.results, types,
                                  feature.names = feature.names,
                                  rank.metric = DMEA.rank.var, p = p, FDR = FDR,
                                  num.permutations = num.permutations,
                                  stat.type = stat.type, 
                                  min.per.set = min.per.set,
                                  scatter.plots = scatter.plots,
                                  scatter.plot.type = scatter.plot.type,
                                  n.dot.sets = n.dot.sets
        ) 
      } else if (DMEA.type == "gene_corr") {
        DMEA.results <- 
          panSEA::mDMEA_gene_corr(drug.sensitivity, gmt.drugs, expression, 
                                  DEGs$all.results, types,
                                  feature.names = feature.names,
                                  rank.metric = DMEA.rank.var, p = p, FDR = FDR,
                                  num.permutations = num.permutations,
                                  stat.type = stat.type, 
                                  min.per.set = min.per.set,
                                  scatter.plots = scatter.plots,
                                  scatter.plot.type = scatter.plot.type,
                                  n.dot.sets = n.dot.sets
          ) 
      }

      # compile inputs & outputs for network graph
      inputs <- list()
      for (i in 1:length(types)) {
        inputs[[types[i]]] <- DMEA.results$all.results[[types[i]]]$corr.result
      }

      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- DMEA.results$all.results[[types[i]]]$result
      }

      DMEA.network <- panSEA::netSEA(
        inputs, outputs, rep("Drug", length(inputs)),
        DMEA.rank.var, p, FDR, n.network.sets, scale
      )
    } else {
      DMEA.results <- NA
      DMEA.network <- NA
    }
  } else if (length(group.names) == 1) {
    # make sure feature names are in column 1
    for (i in 1:length(types)) {
      data.list[[types[i]]] <-
        data.list[[i]][, c(feature.names[i], 
                           colnames(data.list[[i]])
                           [colnames(data.list[[i]]) != feature.names[i]])]
    }
    
    ssGSEA.results <- list()
    ssGSEA.network <- list()
    DMEA.results <- list()
    DMEA.network <- list()
    if (length(group.samples[[1]]) > 1) {
      group.samples <- group.samples[[1]]
    }
    for (j in group.samples) {
      ## extract omics data for each sample
      temp.data <- list()
      for (i in 1:length(types)) {
        temp.sample <- colnames(data.list[[i]])[j]
        temp.data[[types[i]]] <-
          data.list[[i]][, c(feature.names[i], temp.sample)]
      }

      ## ssGSEA & network graph
      if (GSEA & !is.null(gmt.features)) {
        ssGSEA.results[[colnames(data.list[[1]])[j]]] <-
          panSEA::mGSEA(temp.data, gmt.features, types, feature.names,
                        GSEA.rank.var, p = p, FDR = FDR, 
                        num.permutations = num.permutations, 
                        stat.type = stat.type, min.per.set = min.per.set,
                        n.dot.sets = n.dot.sets
          )

        # compile inputs & outputs for network graph
        outputs <- list()
        for (i in 1:length(types)) {
          outputs[[types[i]]] <- 
            ssGSEA.results[[colnames(data.list[[1]])[j]]]$all.results[[types[i]]]$result
        }

        ssGSEA.network[[colnames(data.list[[1]])[j]]] <-
          panSEA::netSEA(
            temp.data, outputs, feature.names,
            GSEA.rank.var, p, FDR, n.network.sets, scale
          )
      } else {
        ssGSEA.results[[colnames(data.list[[1]])[j]]] <- NA
        ssGSEA.network[[colnames(data.list[[1]])[j]]] <- NA
      }

      ## DMEA & network graph
      if (DMEA & !is.null(drug.sensitivity) & !is.null(expression) &
        !is.null(gmt.drugs)) {
        if (DMEA.type == "WV") {
          DMEA.results[[colnames(data.list[[1]])[j]]] <-
            panSEA::mDMEA(drug.sensitivity, gmt.drugs, expression, temp.data,
                          types, feature.names = feature.names, 
                          rank.metric = DMEA.rank.var, 
                          weight.values = GSEA.rank.var, p = p, FDR = FDR, 
                          num.permutations = num.permutations, 
                          stat.type = stat.type, 
                          min.per.set = min.per.set, 
                          scatter.plots = scatter.plots,
                          scatter.plot.type = scatter.plot.type, 
                          n.dot.sets = n.dot.sets
            ) 
        } else if (DMEA.type == "cell_corr") {
          DMEA.results[[colnames(data.list[[1]])[j]]] <-
            panSEA::mDMEA_cell_corr(drug.sensitivity, gmt.drugs, expression, 
                                    temp.data, types, 
                                    feature.names = feature.names, 
                                    rank.metric = DMEA.rank.var, p = p, 
                                    FDR = FDR, 
                                    num.permutations = num.permutations, 
                                    stat.type = stat.type, 
                                    min.per.set = min.per.set, 
                                    scatter.plots = scatter.plots,
                                    scatter.plot.type = scatter.plot.type, 
                                    n.dot.sets = n.dot.sets
            )
        } else if (DMEA.type == "gene_corr") {
          DMEA.results[[colnames(data.list[[1]])[j]]] <-
            panSEA::mDMEA_gene_corr(drug.sensitivity, gmt.drugs, expression, 
                                    temp.data, types, 
                                    feature.names = feature.names, 
                                    rank.metric = DMEA.rank.var, p = p, 
                                    FDR = FDR, 
                                    num.permutations = num.permutations, 
                                    stat.type = stat.type, 
                                    min.per.set = min.per.set, 
                                    scatter.plots = scatter.plots,
                                    scatter.plot.type = scatter.plot.type, 
                                    n.dot.sets = n.dot.sets
            )
        }

        # compile inputs & outputs for network graph
        inputs <- list()
        for (i in 1:length(types)) {
          inputs[[types[i]]] <- 
            DMEA.results[[colnames(data.list[[1]])[j]]]$all.results[[types[i]]]$corr.result
        }

        outputs <- list()
        for (i in 1:length(types)) {
          outputs[[types[i]]] <- 
            DMEA.results[[colnames(data.list[[1]])[j]]]$all.results[[types[i]]]$result
        }

        DMEA.network[[colnames(data.list[[1]])[j]]] <-
          panSEA::netSEA(
            inputs, outputs, rep("Drug", length(inputs)),
            DMEA.rank.var, p, FDR, n.network.sets, scale
          )
      } else {
        DMEA.results[[colnames(data.list[[1]])[j]]] <- NA
        DMEA.network[[colnames(data.list[[1]])[j]]] <- NA
      }
    }
  }

  return(list(
    mDEG.results = DEGs,
    mGSEA.results = ssGSEA.results, mDMEA.results = DMEA.results,
    mGSEA.network = ssGSEA.network, mDMEA.network = DMEA.network
  ))
}
