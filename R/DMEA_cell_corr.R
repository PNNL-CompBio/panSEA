DMEA_cell_corr <- function(drug.sensitivity, gmt = NULL, expression, weights,
                           value = "AUC", sample.names = colnames(expression)[1], 
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
  # reformat expression data frame with cell lines as column names
  rownames(expression) <- expression[ , sample.names]
  expr <- as.data.frame(t(expression))
  expr <- expr[rownames(expr) %in% weights[ , gene.names],]
  expr[ , gene.names] <- rownames(expr)
  
  # merge expression data frame with input weights
  expr.weights <- merge(weights[ , c(gene.names, weight.values)], 
                        expr, by = gene.names)
  
  # for each cell line: run correlation between expression & input weights
  # data points are genes
  expr.weights.corr <- DMEA::rank_corr(expr.weights, variable = sample.names, 
                               plots = FALSE)$result
  
  # merge expr.weights.corr with drug sensitivity
  drug.expr.corr <- merge(expr.weights.corr[ , c(sample.names, rank.metric)], 
                          drug.sensitivity, by = sample.names)
  
  # for each drug: run correlation between sensitivity & drug.expr.corr
  # data points are cell lines
  drug.corr <- DMEA::rank_corr(drug.expr.corr, variable = drug, value = value,
                               type = scatter.plot.type, 
                               min.per.corr = min.per.corr,
                               plots = scatter.plots, FDR = FDR.scatter.plots, 
                               xlab = xlab, ylab = ylab, 
                               position.x = position.x, position.y = position.y,
                               se = se)
  
  # run drugSEA
  if (ties) {
    DMEA.results <- panSEA::drugSEA_ties(
      drug.corr$result, gmt, drug, rank.metric, set.type, FDR = FDR,
      num.permutations = num.permutations, stat.type = stat.type,
      min.per.set = min.per.set, sep = sep, exclusions = exclusions,
      descriptions = descriptions
    )
  } else {
    DMEA.results <- DMEA::drugSEA(
    drug.corr$result, gmt, drug, rank.metric, set.type, FDR = FDR,
    num.permutations = num.permutations, stat.type = stat.type,
    min.per.set = min.per.set, sep = sep, exclusions = exclusions,
    descriptions = descriptions
  )
  }

  return(list(
    cell.corr.result = expr.weights.corr,
    unused.weights = weights[!(weights[ , gene.names] %in% expr[ , gene.names]), ],
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
