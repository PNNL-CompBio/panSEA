DMEA_gene_corr <- function(gmt = NULL, expression, weights,
                           value = "Perturbation", sample.names = colnames(expression)[1], 
                           gene.names = colnames(weights)[1], 
                           weight.values = colnames(weights)[2],
                           rank.metric = "Pearson.est", FDR = 0.25,
                           num.permutations = 1000, stat.type = "Weighted",
                           drug.info = NULL, drug = "Drug", set.type = "moa",
                           min.per.set = 6, sep = "[|]",
                           exclusions = c("-666", "NA", "na", "NaN", "NULL"),
                           descriptions = NULL, min.per.corr = 3, 
                           scatter.plots = TRUE,
                           scatter.plot.type = "pearson", 
                           FDR.scatter.plots = 0.05,
                           xlab = "Expression Correlation Estimate", 
                           ylab = value,
                           position.x = "min", position.y = "min", se = TRUE, ties = FALSE) {
  # merge expression data frame with input weights
  expr.weights <- merge(weights[ , c(gene.names, weight.values)], 
                        expression, by = gene.names)
  
  # for each drug line: run correlation between expression & input weights
  # data points are genes
  drug.corr <- DMEA::rank_corr(expr.weights, variable = sample.names,
                                       value = value, type = scatter.plot.type, 
                                       min.per.corr = min.per.corr,
                                       plots = scatter.plots, 
                                       FDR = FDR.scatter.plots, 
                                       xlab = xlab, ylab = ylab, 
                                       position.x = position.x, 
                                       position.y = position.y, se = se)
  
  # run drugSEA
  DMEA.results <- panSEA::drugSEA_ties(
    drug.corr$result, gmt, drug, rank.metric, set.type, FDR = FDR,
    num.permutations = num.permutations, stat.type = stat.type,
    min.per.set = min.per.set, sep = sep, exclusions = exclusions,
    descriptions = descriptions, ties = ties
  )

  return(list(
    corr.result = drug.corr$result,
    corr.scatter.plots = drug.corr$scatter.plots,
    gmt = gmt,
    result = DMEA.results$result,
    mtn.plots = DMEA.results$mtn.plots,
    volcano.plot = DMEA.results$volcano.plot,
    removed.sets = DMEA.results$removed.sets,
    unannotated.drugs = DMEA.results$unannotated.drugs
  ))
}
