\name{mDMEA}

\alias{mDMEA}

\title{Multi-omic Drug Mechanism Enrichment Analysis
}

\description{For each input omics type: ranks samples using Weighted Voting,
correlates Weighted Voting scores to drug sensitivity,
and then runs enrichment analysis using correlation estimates to
determine the enrichment of drug mechanisms-of-action.
}

\usage{
mDMEA(drug.sensitivity="PRISM", gmt="PRISM", 
  expression=as.list(rep("adherent CCLE", length(types))), weights, types, 
  value="AUC", sample.names=colnames(drug.sensitivity)[1], 
  feature.names=rep("Gene", length(types)), 
  weight.values=rep("Log2FC", length(types)), 
  rank.metric=rep("Pearson.est", length(types)), p=0.05, FDR=0.25, 
  num.permutations=1000, stat.type="Weighted", drug.info=NULL, drug="Drug", 
  set.type="moa", min.per.set=6, sep="[|]", 
  exclusions=c("-666", "NA", "na", "NaN", "NULL"), 
  descriptions=NULL, min.per.corr=3, scatter.plots=TRUE, 
  scatter.plot.type="pearson", FDR.scatter.plots=0.05, 
  xlab="Weighted Voting Score", ylab=value, position.x="min", 
  position.y="min", se=TRUE, n.dot.sets=10)
}

\arguments{
  \item{drug.sensitivity}{Data frame containing drug sensitivity scores
  (e.g., AUC) with drug names as column names, except for one column which 
  contains sample names. If "PRISM", drug AUC sensitivity scores for cancer 
  cell lines from the PRISM drug screen will be used. Default: "PRISM".
  }

  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from
  GSA.read.gmt function from GSA package). If gmt is not provided, drug.info
  parameter must be provided. If "PRISM", drug mechanisms of action 
  annotations from the PRISM drug screen will be used. Default: "PRISM".
  }

  \item{expression}{List of data frames containing normalized expression with 
  feature names (e.g., gene symbols) as column names, except for one column 
  which contains sample names. If "adherent CCLE", normalized RNA-seq data for 
  327 adherent cancer cell lines from CCLE version 19Q4 will be used. Default: 
  as.list(rep("adherent CCLE", length(types))).
  }

  \item{weights}{List of data frames with feature.names and weight.values (e.g.,
  log2(fold-change) values).
  }

  \item{types}{Vector of character descriptions of each omics type (e.g., 
  c("Transcriptomics", "Proteomics")).
  }

  \item{value}{Description of the values starting at column 3. Do not include
  spaces! Default: "AUC".
  }

  \item{sample.names}{Column name in drug sensitivity and expression data
  frames containing the sample names. Default: colnames(drug.sensitivity)[1].
  }

  \item{feature.names}{Vector of column names in weights data frames containing 
  gene symbols. Default: rep("Gene", length(types)).
  }

  \item{weight.values}{Column name in weights data frame containing weight
  values (e.g., log2(fold-change) values). Default: rep("Log2FC", length(types)).
  }

  \item{rank.metric}{Vector of estimates for ranking drugs for enrichment 
  analysis. Possible values include "Pearson.est", "Spearman.est", "Slope", 
  "Rank.slope", or another column name in the output produced by the rank_corr 
  function. Default: rep("Pearson.est", length(types)).
  }

  \item{p}{P-value cut-off for significant enrichments. Default: 0.05.
  }
  
  \item{FDR}{False discovery rate cut-off value for significant enrichments.  
  Default: 0.25.
  }

  \item{num.permutations}{Number of permutations to perform for enrichment
  analysis. Default: 1000.
  }

  \item{stat.type}{Statistic type for enrichment analysis. Must be either
  "Weighted" (score weight = 1) or "Classic" (score weight = 0). Default:
  "Weighted".
  }

  \item{drug.info}{Data frame with set membership information for each element.
  Optional unless gmt is not provided. Default: NULL.
  }

  \item{drug}{Column name in drug.info containing the drug names. Optional
  unless gmt is not provided. Default: "Drug".
  }

  \item{set.type}{Column name in drug.info containing the set membership
  information for each drug. Optional unless gmt is not provided. Default:
  "moa".
  }

  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.
  }

  \item{sep}{Character(s) to use for splitting drug set names. Optional unless
  gmt is not provided. Default: "[|]".
  }

  \item{exclusions}{A vector of drug set names to be excluded from the generated
  gmt. Optional unless gmt is not provided. Default: 
  c("-666", "NA", "na", "NaN", "NULL").
  }

  \item{descriptions}{Optional: column name in drug.info containing
  descriptions for each drug set. If not provided and gmt is also not provided,
  set.type will be used for the descriptions instead. Default: NULL.
  }

  \item{min.per.corr}{Minimum number of samples required to run correlations.
  Default: 3.
  }

  \item{scatter.plots}{Logical operator describing whether or not to produce
  scatter plots for correlations. Default: TRUE.
  }

  \item{scatter.plot.type}{Correlation type for plots. Must be specified as
  "pearson" or "spearman". Default: "pearson".
  }

  \item{FDR.scatter.plots}{False discovery rate cut-off value to produce
  scatter plots for significant correlations. Default: 0.05.
  }

  \item{xlab}{Text describing rank metric for x-axis label of scatter plots.
  Default: "Weighted Voting Score".
  }

  \item{ylab}{Text describing variable values for y-axis label of scatter
  plots. Default: value.
  }

  \item{position.x}{X-axis position of correlation estimate on scatter plots.
  Must be "min", "mid", "max", or numeric. Default: "mid".
  }

  \item{position.y}{Y-axis position of correlation estimate on scatter plots.
  Must be "min", "mid", "max", or numeric. Default: "max".
  }

  \item{se}{Logical operator describing whether or not to include the 95
  percent confidence interval limits of the regression line on the scatter
  plots. Default: TRUE.
  }
  
  \item{n.dot.sets}{Number of sets from enrichment analyses to be represented 
  in the dot plot. The top significant sets will be selected based on 
  abs(mean_NES) and p & FDR cutoff values. Default: 10.
  }
}

\value{
  \item{compiled.results}{Outputs from compile_mDMEA function.}

  \item{all.results}{List of outputs from DMEA function for each input type.}
}

\references{
Garana, B.B., Joly, J.H., Delfarah, A. et al. Drug mechanism enrichment analysis
improves prioritization of therapeutics for repurposing. BMC Bioinformatics 24,
215 (2023). https://doi.org/10.1186/s12859-023-05343-8
}

\author{Belinda B. Garana
}

\section{Warning}{The column name for your samples must be the same in your
drug sensitivity data frame as in your expression data frame. You can set the
name of the column containing your samples using the parameter sample.names.

If no gmt is provided, then drug.info with set.type annotations must be
provided. Either the gmt or set.type column must provide annotations for 2+
drug sets.

Each drug must have only one numeric, nonzero rank.
}


\note{Output venn.diagram in compiled.results will only display results for 
the first 4 data frames in the input data.list.
}

\seealso{
  \link[mDMEA]{compilDMEA}
}

\examples{
## Step 1: prepare drug sensitivity data frame
# create list of sample names
Sample_ID <- seq(from = 1, to = 21)

drug.sensitivity <- as.data.frame(Sample_ID)

# create list of drug names
Drug <- paste0("Drug_", seq(from = 1, to = 24))

# give each drug values representative of AUC sensitivity scores
for(i in 1:length(Drug)){
  drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID),
                                        mean = 0.83, sd = 0.166)
}

## Step 2: prepare drug info data frame
# alternatively, a gmt object could be provided
info <- as.data.frame(Drug)

# moa is our default name for the column containing drug set annotations
info$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 4)]), 6)

## Step 3: prepare gene weight data frames
# create list of gene symbols
Gene <- paste0("Gene_", seq(from = 1, to = 50))

# by default, gene symbols are found in
# the first column of your weights data frame
weights <- as.data.frame(Gene)
weights2 <- as.data.frame(Gene)

# give each gene a weight
# by default, gene weight values are found in
# the second column of your weights data frame
weights$Log2FC <- rnorm(length(Gene))
weights2$Log2FC <- rnorm(length(Gene))

## Step 4: prepare expression data frames
# by default, column 1 of your expression data frame is
# the column name from which sample names are gathered and
# the column containing sample names in your drug sensitivity
# data frame should have the same name
expr <- as.data.frame(Sample_ID)
expr2 <- as.data.frame(Sample_ID)

# give each gene values representative of normalized expression
# each gene is represented by a column in your expression data frame
for(i in 1:length(Gene)){
  expr[,c(Gene[i])] <- rnorm(length(Sample_ID), sd = 0.5)
  expr2[,c(Gene[i])] <- rnorm(length(Sample_ID), sd = 0.5)
}

## Step 5: perform DMEA and store results
mDMEA.test <- panSEA::mDMEA(drug.sensitivity, gmt = NULL, 
  expression = list(expr, expr2), 
  types = c("Transcriptomics", "Proteomics"), 
  weights = list(weights, weights2), drug.info = info)
}
