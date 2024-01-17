\name{panSEA}

\alias{panSEA}

\title{Multi-omic Enrichment Analysis
}

\description{For each input omics type: runs differential expression analysis 
if appropriate followed by enrichment analyses of feature sets and drug 
mechanisms of action. Also generates network graphs of top significantly 
enriched feature and drug sets.
}

\usage{
panSEA(data.list, types, feature.names = rep("Gene", length(types)), 
  GSEA.rank.var = rep("Log2FC", length(types)), 
  DMEA.rank.var = rep("Pearson.est", length(types)),
  group.names=c("Diseased", "Healthy"), 
  group.samples=list(2:(0.5 * (ncol(data.list[[1]]) + 1)),
  (0.5 * (ncol(data.list[[1]]) + 1) + 1):ncol(data.list[[1]])), 
  gmt.features=as.list(rep("msigdb_Homo sapiens_C2_CP:KEGG", length(types))), 
  gmt.drugs="PRISM", GSEA=TRUE, DMEA=TRUE, DMEA.type="WV", p=0.05, FDR=0.25, 
  FDR.features=0.05, min.per.set=6, drug.sensitivity="PRISM", 
  expression=as.list(rep("adherent CCLE", length(types))), 
  n.network.sets=2*length(types), n.dot.sets=10, scale=5)
}

\arguments{
  \item{data.list}{List of data frames where each data frame has one column 
  containing feature names (e.g., gene symbols) and the rest of its columns 
  having sample names as column names.
  }
  
  \item{types}{Vector of character descriptions of each omics type (e.g., 
  c("Transcriptomics", "Proteomics")).
  }
  
  \item{feature.names}{Vector of column names containing feature names (e.g., 
  gene symbols) in each corresponding data frame in the input data.list. 
  Default: rep("Gene", length(inputs)).
  }
  
  \item{GSEA.rank.var}{Vector of column names containing rank variable for 
  gene set enrichment analysis for each corresponding input data frame. If 
  group.names has length 2, column name must be an output of differential 
  expression analysis. If group.names has length 1, column name must match a 
  column name in the corresponding data frame in the input data.list.
  Default: rep("Log2FC", length(types)).
  }
  
  \item{DMEA.rank.var}{Vector of column names containing rank variable for 
  drug mechanism enrichment analysis for each corresponding data frame in the 
  input data.list. Column names must match a column name from the result of 
  DMEA::rank_corr. Default: rep("Pearson.est", length(types)).
  }
  
  \item{group.names}{Vector of names for 2 sample groups. 
  Default: c("Diseased", "Healthy").
  }
  
  \item{group.samples}{List of column numbers or names corresponding to the 
  group names. For example, if columns 2:4 correspond to the first group (e.g., 
  "Diseased") and columns 5:7 correspond to the second group (e.g., "Healthy"), 
  then group.samples would be list(2:4, 5:7). The order of samples must be the 
  same across all data frames in the input data.list. Default: 
  list(2:(0.5 * (ncol(data.list[[1]]) + 1)),
  (0.5 * (ncol(data.list[[1]]) + 1) + 1):ncol(data.list[[1]])).
  }
  
  \item{gmt.features}{List of gmt objects containing sets of features for gene 
  set enrichment analyses (e.g., output from GSA.read.gmt function from GSA 
  package)  or character strings with input parameters for msigdbr::msigdbr 
  separated by "_". Default: 
  rep("msigdb_Homo sapiens_C2_CP:KEGG", length(types)).
  }
  
  \item{gmt.drugs}{Gmt object containing drug sets (e.g., output from
  GSA.read.gmt function from GSA package). If "PRISM", drug mechanisms of 
  action annotations from the PRISM drug screen will be used. Default: "PRISM".
  }
  
  \item{GSEA}{Logical operator describing whether or not to run gene set 
  enrichment analysis. Default: TRUE.
  }
  
  \item{DMEA}{Logical operator describing whether or not to run drug mechanism 
  enrichment analysis. Default: TRUE.
  }
  
  \item{DMEA.type}{Character string describing type of drug rank list to use for 
  drug mechanism enrichment analysis. Must be either "WV" for drugs to be 
  ranked by correlations with weighted voting expression scores, "cell_corr" for
  correlations with expression values, or "gene_corr" for correlations with 
  drug expression perturbation values. Default: "WV".
  }
  
  \item{p}{P-value cut-off for significant enrichments. Default: 0.05.
  }

  \item{FDR}{False discovery rate cut-off value for significant enrichments. 
  Default: 0.25.
  }
  
  \item{FDR.features}{False discovery rate cut-off value for significant 
  differential expression. Default: 0.05.
  }
  
  \item{num.permutations}{Number of permutations to perform for enrichment
  analysis. Default: 1000.
  }

  \item{stat.type}{Statistic type for enrichment analysis. Must be either
  "Weighted" (score weight = 1) or "Classic" (score weight = 0). Default:
  "Weighted".
  }
  
  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.
  }
  
  \item{scatter.plots}{Logical operator describing whether or not to produce
  scatter plots for DMEA correlations. Default: TRUE.
  }

  \item{scatter.plot.type}{Correlation type for DMEA scatter plots. Must be 
  specified as "pearson" or "spearman". Default: "pearson".
  }
  
  \item{drug.sensitivity}{Data frame containing drug sensitivity scores
  (e.g., AUC) with drug names as column names, except for one column which 
  contains sample names. If "PRISM", drug AUC sensitivity scores for cancer 
  cell lines from the PRISM drug screen will be used. Default: "PRISM".
  }

  \item{expression}{List of data frames containing normalized expression with 
  feature names (e.g., gene symbols) as column names, except for one column 
  which contains sample names. If "adherent CCLE", normalized RNA-seq data for 
  327 adherent cancer cell lines from CCLE version 19Q4 will be used. 
  Default: as.list(rep("adherent CCLE", length(types))).
  }
  
  \item{n.network.sets}{Number of sets from enrichment analyses to be 
  represented in the network graphs. The top significant sets will be selected 
  based on abs(NES) and p & FDR cutoff values. Default: 2*length(inputs).
  }
  
  \item{n.dot.sets}{Number of sets from enrichment analyses or features from 
  differential expression analyses to be represented in the dot plots. The top 
  significant results will be selected based on abs(NES) for sets and 
  abs(Log2FC) for features as well as p & FDR cutoff values. Default: 10.
  }
  
  \item{scale}{Numeric value to multiply by rank.var value to size vertices on 
  network graphs for visibility. Default: 5.
  }
}

\value{
  \item{mDEG.results}{Outputs from mDEG function. NA if length(group.names) == 
  1.}
  
  \item{mGSEA.results}{Outputs from mGSEA function. If 
  length(group.names) == 1, list of outputs corresponds to each sample (i.e., 
  column name in each data frame in the input data.list)}
  
  \item{mDMEA.results}{Outputs from mDMEA function. If 
  length(group.names) == 1, list of outputs corresponds to each sample (i.e., 
  column name in each data frame in the input data.list)}
  
  \item{mGSEA.network}{Outputs from netSEA function containing network graphs 
  of leading edge features from top significantly enriched sets from mGSEA 
  analysis. If  length(group.names) == 1, list of outputs corresponds to each 
  sample (i.e., column name in each data frame in the input data.list)}
  
  \item{mDMEA.network}{Outputs from netSEA function containing network graphs 
  of leading edge drugs from top significantly enriched sets from mDMEA 
  analysis. If length(group.names) == 1, list of outputs corresponds to each 
  sample (i.e., column name in each data frame in the input data.list)}
}

\references{
Garana, B.B., Joly, J.H., Delfarah, A. et al. Drug mechanism enrichment analysis
improves prioritization of therapeutics for repurposing. BMC Bioinformatics 24,
215 (2023). https://doi.org/10.1186/s12859-023-05343-8
}

\author{Belinda B. Garana
}

\note{To save static network graphs (e.g., mGSEA.network$static and 
mDMEA.network$static outputs), generate plot and then use dev.print. 

To save 
interactive network graphs (e.g., mGSEA.network$interactive and 
mDMEA.network$interactive outputs) as HTML file, use visNetwork::visSave.

For DMEA.type input parameter: "WV" is recommended for an input gene signature 
(e.g., representing Diseased vs. Healthy phenotypes) when you are seeking toxic
drugs, "cell_corr" is recommended for an input expression data frame (e.g., 
expression data for individual samples representing one phenotype) when you are 
seeking toxic drugs, and "gene_corr" is recommended for an input gene signature 
(e.g., representing Diseased vs. Healthy phenotypes) when you are seeking drugs 
with similar/opposite effects on gene expression.
}

\section{Warning}{The column name for your samples must be the same in your
drug sensitivity data frame as in your expression data frame.
}

\seealso{
  \link[panSEA]{mDEG}
  \link[panSEA]{mGSEA}
  \link[panSEA]{mDMEA}
  \link[panSEA]{mDMEA_cell_corr}
  \link[panSEA]{mDMEA_gene_corr}
  \link[panSEA]{netSEA}
}

\examples{
### example of 1 group: signatures of sensitivity to RAF or EGFR inhibitors
## Step 1: get inputs & format
examples <- c("GSE66539_SKCM_sensitive_vs_resistant_to_vemurafenib",
"GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib")
types <- c("GSE66539", "GSE31625") # short names for inputs for plots

inputs <- list()
for (i in 1:length(examples)) {
inputs[[examples[i]]] <- read.csv(
paste0("https://raw.githubusercontent.com/BelindaBGarana/", 
"DMEA/shiny-app/Examples/Gene_signature/", examples[i], 
"/Filtered_gene_signature_no_duplicates.csv"))
}

## Step 2: run multi-omic enrichment analyses
panSEA.test <- panSEA::panSEA(inputs, types, 
group.names = "Sensitive_vs_resistant_to_inhibitor", group.samples = list(2))
}
