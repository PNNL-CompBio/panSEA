panSEA <- function(data, types, group.names="Diseased", 
                   group.samples=1:ncol(data[[1]]), 
                   gmt.features="msigdb_Homo sapiens_C2_CP:KEGG", 
                   gmt.drugs="PRISM", FDR=0.25, min.per.set=6, 
                   drug.sensitivity="PRISM", expression="adherent CCLE",
                   n.top.sets = length(types)){
  #### Step 1. Check if formats are correct ####
  # make sure that there are as many types as gmt.features inputs
  if (is.character(gmt.features)) {
    # get gmt.features if not provided
    if (grepl("msigdb", gmt.features, ignore.case = TRUE)) {
      gmt.info <- stringr::str_split(gmt, "_")[[1]]
      if (length(gmt.info) > 1) {
        if (length(gmt.info) == 2) {
          msigdb.info <- msigdbr::msigdbr(gmt.info[2])
        } else if (length(gmt.info) == 3) {
          msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3])
        } else {
          msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3], gmt.info[4])
        }
        gmt.features <- DMEA::as_gmt(msigdb.info, "gene_symbol", "gs_name", 
                            "gs_description")
      }
    }
    gmt.features <- as.list(rep(gmt.features, length(types))) 
  } else if (length(types) != length(gmt.features)) {
      stop("Length of types vector must match that of gmt.features list")
  }

  if (length(conditions) > 2 | length(conditions) < 1) {
    stop("Only 1 or 2 conditions are allowed")
  } else if (length(conditions) == 2) {
    #### Step 2. Differential expression analysis ####
    deg <- mDEG(data, types, group.names, group.samples)
      
    #### Step 3. Enrichment analyses ####
    ### ssGSEA & network graph
    if (!is.null(gmt.features)) {
      ssGSEA.results <- mGSEA(deg$DEGs, gmt.features, FDR = FDR, 
                              min.per.set = min.per.set)
      
      # compile inputs & results for network graph
      inputs <- data.table::rbindlist(deg$DEGs, use.names = TRUE, 
                                      idcol = "type")
      
      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- ssGSEA.results[[types[i]]]$result
      }
      outputs <- data.table::rbindlist(outputs, use.names = TRUE, 
                                       idcol = "type")
      
      ssGSEA.network <- netEA(inputs, outputs, n.top.sets)
    }
    
    ### DMEA & network graph
    if (!is.null(drug.sensitivity) & !is.null(expression) 
        & !is.null(gmt.drugs)) {
      DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, 
                            expression, deg, FDR = FDR)
      
      # compile inputs & outputs for network graph
      inputs <- list()
      for (i in 1:length(types)) {
        inputs[[types[i]]] <- DMEA.results[[types[i]]]$corr.result
      }
      inputs <- data.table::rbindlist(inputs, use.names = TRUE, 
                                      idcol = "type")
      
      outputs <- list()
      for (i in 1:length(types)) {
        outputs[[types[i]]] <- DMEA.results[[types[i]]]$result
      }
      outputs <- data.table::rbindlist(outputs, use.names = TRUE, 
                                       idcol = "type")
      
      DMEA.network <- netEA(inputs, outputs, n.top.sets)
    } 
  } else if (length(conditions) == 1) {
    for (j in 2:ncol(data)) {
      #### Step 2. Extract omics data for each sample ####
      temp.data <- list()
      for (i in 1:length(types)) {
        temp.data[[types[i]]] <- data[[types[i]]][ , c(1,j)]
      }
      
      #### Step 3. Enrichment analyses ####
      ### ssGSEA & network graph
      if (!is.null(gmt.features)) {
        ssGSEA.results <- mGSEA(temp.data, gmt.features, FDR = FDR, 
                                min.per.set = min.per.set)
        
        # compile inputs & outputs for network graph
        inputs <- data.table::rbindlist(temp.data, use.names = TRUE, 
                                        idcol = "type")
        
        outputs <- list()
        for (i in 1:length(types)) {
          outputs[[types[i]]] <- ssGSEA.results[[types[i]]]$result
        }
        outputs <- data.table::rbindlist(outputs, use.names = TRUE, 
                                         idcol = "type")
        
        ssGSEA.network <- netEA(inputs, outputs, n.top.sets)
      }
      
      ### DMEA ###
      if (!is.null(drug.sensitivity) & !is.null(expression) 
          & !is.null(gmt.drugs)) {
        DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, expression, 
                              temp.data, FDR = FDR, min.per.set = min.per.set)
        
        # compile inputs & outputs for network graph
        inputs <- list()
        for (i in 1:length(types)) {
          inputs[[types[i]]] <- DMEA.results[[types[i]]]$corr.result
        }
        inputs <- data.table::rbindlist(inputs, use.names = TRUE, 
                                        idcol = "type")
        
        outputs <- list()
        for (i in 1:length(types)) {
          outputs[[types[i]]] <- DMEA.results[[types[i]]]$result
        }
        outputs <- data.table::rbindlist(outputs, use.names = TRUE, 
                                         idcol = "type")
        
        DMEA.network <- netEA(inputs, outputs, n.top.sets)
      }
    }
  }

  return(list(DEGs = deg$DEGs, Log2Transformed = deg$Log2Transformed, 
              ssGSEA.results = ssGSEA.results, DMEA.results = DMEA.results, 
              ssGSEA.network = ssGSEA.network, DMEA.network = DMEA.network))
}
