mGSEA <- function(data, gmt, types, direction.adjust=NULL, FDR=0.25, 
                  num.permutations=1000, stat.type="Weighted", min.per.set=6){
  #### Step 1. Perform ssGSEA on each omics type ####
  ssGSEA.list <- list()
  for(i in 1:length(types)) {
    message(paste("Running ssGSEA using", types[i], "data"))
    
    ssGSEA.list[[types[i]]] <- ssGSEA(data[[i]], 
                                    gmt[[i]], direction.adjust, FDR,
                                    num.permutations, stat.type, min.per.set,
                                    sep, exclusions, descriptions)
  }

  #### Step 2. Compile ssGSEA results across omics types ####
  compiled.GSEA <- compilGSEA(ssGSEA.list)
  
  return(list(compiled.results = compiled.GSEA, 
              all.results = ssGSEA.list
  ))
}
