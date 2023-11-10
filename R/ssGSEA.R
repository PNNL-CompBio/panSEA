ssGSEA <- function(data, gmt=NULL, direction.adjust=NULL, FDR=0.25, 
                   num.permutations=1000, stat.type="Weighted", min.per.set=6, 
                   sep = ", ", 
                   exclusions = c("-666", "NA", "na", "NaN", "NULL"), 
                   descriptions=NULL){
  
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
