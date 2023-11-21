\name{mDEG}

\alias{mDEG}

\title{Multi-omic Diffential Expression Analysis
}

\description{Run differential expression analysis for each input data type.
}

\usage{
mDEG(data, types, group.names=c("Diseased", "Healthy"), group.samples)
}

\arguments{
  \item{data}{List of data frames where each data frame has feature names 
  (e.g., gene symbols) as row names and sample names as column names.
  }

  \item{types}{Vector of character descriptions of each omics type (e.g., 
  c("Transcriptomics", "Proteomics")).
  }
  
  \item{group.names}{Vector of names for 2 sample groups. 
  Default: c("Diseased", "Healthy").
  }
  
  \item{group.samples}{Vector of column numbers or names corresponding to the 
  group names. For example, if columns 1:3 correspond to the first group (e.g., 
  "Diseased") and columns 4:6 correspond to the second group (e.g., "Healthy"), 
  then group.samples would be c(1:3, 4:6). The order of samples must be the 
  same across all input data frames.
  }
}

\value{
  \item{DEGs}{List of data frames containing differential expression results 
  after using limma::eBayes for each input data frame.}

  \item{Log2Transformed}{Vector of logical operators describing if the 
  corresponding input data was log2-transformed to normalize the distribution.}
}

\author{Belinda B. Garana
}

\seealso{
  \link[mDEG]{compilDMEA}
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
expr <- as.data.frame(Sample_1, Sample_2, Sample_3, 
  Sample_4, Sample_5, Sample_6)
rownames(expr) <- Gene

expr2 <- as.data.frame(Sample_1, Sample_2, Sample_3, 
  Sample_4, Sample_5, Sample_6)
rownames(expr2) <- Gene

## Step 3: perform DMEA and store results
mDEG.test <- mDEG(list(expr, expr2), 
  types = c("Transcriptomics", "Proteomics"), group.samples = c(1:3, 4:6))
}
