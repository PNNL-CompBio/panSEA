\name{mDEG}

\alias{mDEG}

\title{Multi-omic Diffential Expression Analysis
}

\description{Run differential expression analysis using limma::eBayes for each 
input data type. If the input data is not normally distributed, Log2 
transformation is performed prior to analysis.
}

\usage{
mDEG(data.list, types, group.names=c("Diseased", "Healthy"), 
  group.samples=list(2:(0.5 * (ncol(data.list[[1]]) + 1)),
  (0.5 * (ncol(data.list[[1]]) + 1) + 1):ncol(data.list[[1]])), 
  feature.names=rep("Gene", length(types)), 
  p=0.05, FDR.features=0.05, n.dot.features=10)
}

\arguments{
  \item{data.list}{List of data frames where each data frame has one column 
  containing feature names (e.g., gene symbols) and the rest of its columns 
  having sample names as column names.
  }

  \item{types}{Vector of character descriptions of each omics type (e.g., 
  c("Transcriptomics", "Proteomics")).
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
  
  \item{feature.names}{Vector of column names containing feature names (e.g., 
  gene symbols) in each corresponding data frame in the input data.list. 
  Default: rep("Gene", length(inputs)).
  }
  
  \item{p}{P-value cut-off for significant differential expression. Default: 
  0.05.
  }
  
  \item{FDR.features}{False discovery rate cut-off value for significant 
  differential expression. Default: 0.05.
  }

  \item{n.dot.features}{Number of features from differential expression 
  analyses to be represented in the dot plots. The top significant features 
  will be selected based on abs(Log2FC) and p & FDR cutoff values. Default: 10.
  }
}

\value{ 
  \item{compiled.results}{Outputs from compile_mDEG function.}

  \item{all.results}{List of data frames containing differential expression 
  analysis results for each input type. The last column of each output data 
  frame 'Log2Transformed' describes if the corresponding input data was 
  log2-transformed to normalize the distribution}
}

\author{Belinda B. Garana
}

\examples{
## Step 1: prepare expression data frames
# create list of gene symbols
Gene <- paste0("Gene_", seq(from = 1, to = 50))

# prepare expression data for samples
Sample_1 <- rnorm(length(Gene), sd = 0.5)
Sample_2 <- rnorm(length(Gene), sd = 0.5)
Sample_3 <- rnorm(length(Gene), sd = 0.5)
Sample_4 <- rnorm(length(Gene), sd = 0.5)
Sample_5 <- rnorm(length(Gene), sd = 0.5)
Sample_6 <- rnorm(length(Gene), sd = 0.5)

# feature names (e.g., gene symbols) should be row names and 
# column names should be sample names
expr <- data.frame(Gene, Sample_1, Sample_2, Sample_3, 
  Sample_4, Sample_5, Sample_6)

expr2 <- data.frame(Gene, Sample_1, Sample_2, Sample_3, 
  Sample_4, Sample_5, Sample_6)

## Step 3: perform differential expression analysis
# by default, we assume the first half of samples (1-3) are "Diseased"
# and the second half of samples (4-6) are "Healthy"
mDEG.test <- panSEA::mDEG(list(expr, expr2), 
  types = c("Transcriptomics", "Proteomics"))
}
