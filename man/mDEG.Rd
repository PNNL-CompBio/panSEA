\name{mDEG}

\alias{mDEG}

\title{Multi-omic Diffential Expression Analysis
}

\description{Run differential expression analysis using limma::eBayes for each 
input data type. If the input data is not normally distributed, Log2 
transformation is performed prior to analysis.
}

\usage{
mDEG(data.list, factor.info, feature.names = rep("Gene", length(data.list)),
  p=0.05, FDR.features=0.05, n.dot.features=10)
}

\arguments{
  \item{data.list}{List of data frames where each data frame has one column 
  containing feature names (e.g., gene symbols) and the rest of its columns 
  having sample names as column names. List names will be interpreted as 
  omics types (e.g., "global" or "phospho"). Non-numeric column will be 
  assumed to be feature names (e.g., "Gene" or "SUB_SITE").
  }

  \item{factor.info}{Data frame where rownames are sample names found in each 
  data frame of data.list, the first column is a contrast with 2 levels, and 
  optionally there may be a second column to be included in the design matrix 
  (e.g., sex or batch effect).}
  
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
  will be selected based on abs(mean_Log2FC) and p & FDR cutoff values. 
  Default: 10.
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

\note{Output venn.diagram in compiled.results will only display results for 
the first 4 data frames in the input data.list.
}

\seealso{
  \link[panSEA]{compile_mDEG}
  \link[panSEA]{panSEA}
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

# compile expression data into list with descriptive names
data.list <- list("Transcriptomics" = expr, "Proteomics" = expr2)

# prepare factor info for differential expression contrast
factor.info <- data.frame("contrast" = c(rep(0, 3), rep(1, 3)))
rownames(factor.info) <- colnames(expr)[2:ncol(expr)]

## Step 3: perform differential expression analysis
mDEG.test <- panSEA::mDEG(data.list, factor.info)
}
