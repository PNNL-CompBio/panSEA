ssGSEA <- function(data, gmt = "msigdb_Homo sapiens_C2_CP:KEGG",
                   feature.names = colnames(data)[1],
                   rank.var = colnames(data)[2], direction.adjust = NULL,
                   FDR = 0.25, num.permutations = 1000,
                   stat.type = "Weighted", min.per.set = 6, ties = FALSE) {
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

        # extract necessary info into data frame
        msigdb.info <- as.data.frame(msigdb.info[, c(
          "gene_symbol",
          "gs_name",
          "gs_description"
        )])

        gmt <- DMEA::as_gmt(
          msigdb.info, "gene_symbol", "gs_name", min.per.set,
          descriptions = "gs_description"
        )
      }
    }
  }

  # run ssGSEA
  if (ties) {
    results <- panSEA::drugSEA_ties(
      data, gmt, feature.names, rank.var,
      "gs_name", direction.adjust, FDR,
      num.permutations, stat.type, min.per.set
    )
  } else {
    results <- DMEA::drugSEA(
      data, gmt, feature.names, rank.var,
      "gs_name", direction.adjust, FDR,
      num.permutations, stat.type, min.per.set
    ) 
  }

  # change "Drug_set" column names to "Feature_set"
  colnames(results$result)[2] <- "Feature_set"
  colnames(results$result.w.ties)[2] <- "Feature_set"
  colnames(results$removed.sets)[1] <- "Feature_set"
  colnames(results$removed.sets)[2] <- "N_features"

  # remove irrelevant outputs of drugSEA function
  results$replaced.drugs <- NULL

  # change name of output referring to drugs to features
  if (ties) {
    names(results)[11] <- "unannotated.features"
  } else {
    names(results)[6] <- "unannotated.features" 
  }

  return(results)
}
