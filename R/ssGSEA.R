ssGSEA <- function(data, gmt="msigdb_Homo sapiens_C2_CP:KEGG", direction.adjust=NULL, FDR=0.25, 
                   num.permutations=1000, stat.type="Weighted", min.per.set=6){
  # get gmt if not provided
  if (is.character(gmt)) {
    if (grepl("msigdb", gmt, ignore.case = TRUE)) {
      gmt.info <- stringr::str_split(gmt, "_")[[1]]
      if (length(gmt.info) > 1) {
        if (length(gmt.info) == 2) {
          msigdb.info <- msigdbr::msigdbr(gmt.info[2])
        } else if (length(gmt.info) == 3) {
          msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3])
        } else {
          msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3], gmt.info[4])
        }
        
        gmt <- DMEA::as_gmt(msigdb.info, "gene_symbol", "gs_name", 
                            "gs_description")
      }
    }
  }
  
  # run ssGSEA
  results <- DMEA::drugSEA(data, gmt, colnames(data)[1], colnames(data)[2],
                           "gs_name", direction.adjust, FDR,
                           num.permutations, stat.type, min.per.set,
                           sep, exclusions, descriptions)
  
  # change "Drug_set" column names to "Gene_set"
  colnames(results$result)[1] <- "Gene_set"
  colnames(results$removed.sets)[1] <- "Gene_set"
  
  # remove irrelevant outputs of drugSEA function
  results$replaced.drugs <- NULL
  
  return(results)
}
