mGSEA <- function(data.list, gmt = as.list(rep(
                    "msigdb_Homo sapiens_C2_CP:KEGG",
                    length(types)
                  )),
                  types, feature.names = rep("Gene", length(types)),
                  rank.var = rep("Log2FC", length(types)),
                  direction.adjust = NULL, p = 0.05, FDR = 0.25,
                  num.permutations = 1000, stat.type = "Weighted", 
                  min.per.set = 6, n.dot.sets = 10, ties = FALSE) {
  #### Step 1. Perform ssGSEA on each omics type ####
  ssGSEA.list <- list()
  for (i in 1:length(types)) {
    message(paste("Running ssGSEA using", types[i], "data"))

    ssGSEA.list[[types[i]]] <- panSEA::ssGSEA(
      data.list[[i]],
      gmt[[i]], feature.names[i], rank.var[i],
      direction.adjust, FDR,
      num.permutations, stat.type, min.per.set, ties
    )
  }

  #### Step 2. Compile ssGSEA results across omics types ####
  if (length(types) > 1) {
    compiled.GSEA <- panSEA::compile_mGSEA(ssGSEA.list, p, FDR, n.dot.sets)
  } else {
    compiled.GSEA <- NA
  }
  

  return(list(
    compiled.results = compiled.GSEA,
    all.results = ssGSEA.list
  ))
}
