DMEA_gene_corr <- function(drug.sensitivity, gmt = NULL, expression, weights,
                           value = "AUC", sample.names = colnames(expression)[1], 
                           gene.names = colnames(data)[1], 
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
                           position.x = "min", position.y = "min", se = TRUE) {
  # merge expression data frame with input weights
  expr.weights <- merge(weights[ , c(gene.names, weight.values)], 
                        expression, by = gene.names)
  
  # for each drug line: run correlation between expression & input weights
  # data points are genes
  expr.weights.corr <- DMEA::rank_corr(expr.weights, variable = sample.names, 
                               plots = FALSE)$result
  
  # run drugSEA
  results <- DMEA::drugSEA(
    expr.weights.corr, gmt, drug, rank.metric, set.type, FDR = FDR,
    num.permutations = num.permutations, stat.type = stat.type,
    min.per.set = min.per.set, sep = sep, exclusions = exclusions,
    descriptions = descriptions
  )

  return(results)
}
