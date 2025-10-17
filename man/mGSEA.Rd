\name{mGSEA}

\alias{mGSEA}

\title{Multi-omic Gene Set Enrichment Analysis
}

\description{Runs enrichment analysis for each input omics type.
}

\usage{
mGSEA(data.list, gmt=as.list(rep("msigdb_Homo sapiens_C2_CP:KEGG", length(types))), 
  types, feature.names=rep("Gene", length(types)), 
  rank.var=rep("Log2FC", length(types)), direction.adjust=NULL, p=0.05, 
  FDR=0.25, num.permutations=1000, stat.type="Weighted", min.per.set=6, 
  n.dot.sets=10, ties=FALSE)
}

\arguments{
  \item{data.list}{List of data frames where the first column of each data frame 
  contains feature names (e.g., gene symbols) and the second column contains 
  the rank metric (e.g., Log2FC).
  }

  \item{gmt}{List of gmt objects containing sets of features (e.g., output from
  GSA.read.gmt function from GSA package) or character strings with input 
  parameters for msigdbr::msigdbr separated by "_". Default: 
  as.list(rep("msigdb_Homo sapiens_C2_CP:KEGG", length(types))).
  }

  \item{types}{Vector of character descriptions of each omics type (e.g., 
  c("Transcriptomics", "Proteomics")).
  }
  
  \item{feature.names}{Vector of column names containing feature names (e.g., 
  gene symbols) in each corresponding data frame in the input data.list. 
  Default: rep("Gene", length(inputs)).
  }
  
  \item{rank.var}{Vector of column names containing rank variable for 
  enrichment analysis for each corresponding data frame in the input data.list. 
  Default: rep("Log2FC", length(types)).
  }
  
  \item{direction.adjust}{Vector of conditions for set.type to identify which,
  if any, estimates to multiply by -1 (e.g., c(" agonist", " activator")).
  Default: NULL.
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
  
  \item{min.per.set}{The minimum number of drugs required per set. Default: 6.
  }
  
  \item{n.dot.sets}{Number of sets from enrichment analyses to be represented 
  in the dot plot. The top significant sets will be selected based on 
  abs(mean_NES) and p & FDR cutoff values. Default: 10.
  }
  
  \item{ties}{Logical operator describing whether or not to permute order of genes with tied ranks. Default: FALSE.
  }
}

\value{
  \item{compiled.results}{Outputs from compile_mGSEA function.}

  \item{all.results}{List of outputs from ssGSEA function for each input type.}
}

\author{Belinda B. Garana
}

\section{Warning}{Each feature must have only one numeric, nonzero rank.
}

\note{Output venn.diagram in compiled.results will only display results for 
the first 4 data frames in the input data.list.
}

\seealso{
  \link[panSEA]{compile_mGSEA}
  \link[panSEA]{panSEA}
}

\examples{
## Step 1: prepare gene weight data frames
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

# combine input data frames into list
data.list <- list(weights, weights2)

## Step 2: prepare set info
info <- as.data.frame(Gene)

# add column for gene set annotations
info$gs_name <- rep(paste("Set", LETTERS[seq(from = 1, to = 5)]), 10)

# convert data frame into gmt object
gmt <- DMEA::as_gmt(info, element.names = "Gene", set.names = "gs_name")

# create list of gmt objects corresponding to each input in data.list
gmt.list <- list(gmt, gmt)

## Step 3: perform mGSEA and store results
mGSEA.test <- panSEA::mGSEA(data.list, gmt.list, 
                  types = c("Transcriptomics", "Proteomics"))
}
