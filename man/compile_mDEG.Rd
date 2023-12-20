\name{compile_mDEG}

\alias{compile_mDEG}

\title{Compile Multi-omic Differential Expression Analysis Results
}

\description{Compiles enrichment analysis results across different inputs.
}

\usage{
compile_mDEG(DEGs, p=0.05, FDR.features=0.05, n.dot.features=10)
}

\arguments{
  \item{DEGs}{List of data frames containing differential expression analysis 
  results. Each data frame must contain features in column 1 and columns named 
  Log2FC, P.Value, and adj.P.Val.
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
  \item{results}{Data frame of compiled differential expression results.}
  
  \item{mean.results}{Data frame of averaged Log2FC results with Fisher 
  combined p-values adjusted with the Benjamini-Hochberg procedure.}
  
  \item{Log2FC.df}{Data frame of compiled log2(fold-change) values (Log2FC).}
  
  \item{minusLogFDR.df}{Data frame of compiled -log(FDR) values.}
  
  \item{venn.diagram}{ggplot object of venn diagram displaying the number of 
  significantly differentially expressed features across different inputs and 
  percentage of the total number of significantly differentially expressed 
  features belonging to each input.}

  \item{dot.plot}{ggplot2 object of dot plot displaying top significantly 
  differentially expressed features across different inputs.}
  
  \item{corr}{Data frame of correlation results between Log2FC values of 
  different inputs.}
  
  \item{corr.matrix}{Tile plot of correlation results between Log2FC values of 
  different inputs.}
}

\author{Belinda B. Garana
}

\note{Output venn.diagram will only display results for the first 4 data frames 
in the input data.list.
}

\seealso{
  \link[panSEA]{mDEG}
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

## Step 4: compile mDEG results
compile_mDEG.test <- panSEA::compile_mDEG(mDEG.test$all.results)
}
