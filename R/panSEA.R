panSEA <- function(data, types, gmt.genes=NULL, gmt.drugs=NULL, gmt.ptm=NULL, 
                   FDR=0.25, min.per.set=6, drug.sensitivity=NULL, 
                   expression=NULL){
  if (length(conditions) > 2 | length(conditions) < 1) {
    stop("Only 1 or 2 conditions are allowed")
  } else if (length(conditions) == 2) {
    deg <- list()
    deg.ptm <- list()
    #### Step 1. Differential expression analysis ####
    for (i in 1:length(types)) {
      if (types[i] %like% "phospho") {
        deg.ptm[[types[i]]] <- limma::eBayes()
        deg[[types[i]]] <- limma::eBayes()
      } else {
        deg[[types[i]]] <- limma::eBayes()
      }
    }
    
    #### Step 2. Enrichment analyses ####
    ### ssGSEA ###
    if (!is.null(gmt.genes)) {
      ssGSEA.results <- mGSEA(deg, gmt.genes, FDR = FDR, 
                              min.per.set = min.per.set) 
    }
    
    ### KSEA if relevant ###
    if (any(types %like% "phospho") & !is.null(gtm.ptm)) {
      KSEA.results <- mGSEA(deg.ptm, gmt.ptm, FDR=FDR, 
                            min.per.set = min.per.set)
    }
    
    ### DMEA ###
    if (!is.null(drug.sensitivity) & !is.null(expression) 
        & !is.null(gmt.drugs)) {
      DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, 
                            expression, deg, FDR = FDR) 
    } 
    
    #### Step 3. Compile results ####
    compiled <- compilEA(deg, ssGSEA.results, DMEA.results, 
                         deg.ptm, KSEA.results)
    
    #### Step 4. Network graphs ####
    networks <- netEA(compiled)
  } else if (length(conditions) == 1) {
    for (j in 2:ncol(data)) {
      #### Step 1. Extract omics data for each sample ####
      temp.data <- list()
      for (i in 1:length(types)) {
        temp.data[[types[i]]] <- data[[types[i]]][ , c(1,j)]
      }
      
      #### Step 2. Enrichment analyses ####
      ### ssGSEA ###
      if (!is.null(gmt.genes)) {
        ssGSEA.results <- mGSEA(temp.data, gmt.genes, FDR = FDR, 
                                min.per.set = min.per.set) 
      }
      
      ### KSEA if relevant ###
      if (any(types %like% "phospho") & !is.null(gmt.ptm)) {
        KSEA.results <- mGSEA(temp.data, gmt.ptm, FDR=FDR, 
                              min.per.set = min.per.set)
      }
      
      #### Step 3. Compile results ####
      compiled <- compilEA(deg, ssGSEA.results, DMEA.results, 
                           deg.ptm, KSEA.results)
      
      #### Step 3. Network graphs ####
      networks <- netEA(compiled)
    }
    ### DMEA ###
    if (!is.null(drug.sensitivity) & !is.null(expression) 
        & !is.null(gmt.drugs)) {
      DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, expression, 
                            temp.data, FDR = FDR, min.per.set = min.per.set)
    } 
  }
  
  return(list(compiled.results = compiled, network.plots = networks))
}
