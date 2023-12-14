\name{compile_mGSEA}

\alias{compile_mGSEA}

\title{Compile Multi-omic Gene Set Enrichment Analysis Results
}

\description{Compiles enrichment analysis results across different inputs.
}

\usage{
compile_mGSEA(ssGSEA.list, p=0.05, FDR=0.25, n.dot.sets=10)
}

\arguments{
  \item{ssGSEA.list}{List of outputs from different ssGSEA runs.
  }

  \item{p}{P-value cut-off for significant enrichments. Default: 0.05.
  }

  \item{FDR}{False discovery rate cut-off value for significant enrichments. 
  Default: 0.25.
  }
  
  \item{n.dot.sets}{Number of sets from enrichment analyses to be represented 
  in the dot plots. The top significant sets will be selected based on 
  abs(mean_NES) and p & FDR cutoff values. Default: 10.
  }
}

\value{
  \item{results}{Data frame of compiled enrichment results.}
  
  \item{mean.results}{Data frame of averaged enrichment results with Fisher 
  combined p-values adjusted with the Benjamini-Hochberg procedure.}
  
  \item{NES.df}{Data frame of compiled normalized enrichment scores (NES).}
  
  \item{minusLogFDR.df}{Data frame of compiled -log(FDR) values.}
  
  \item{venn.diagram}{ggplot object of venn diagram displaying the number of 
  significantly enriched sets across different inputs and percentage of the 
  total number of significantly enriched sets belonging to each input.}

  \item{dot.plot}{ggplot2 object of dot plot displaying top significantly 
  enriched sets across different inputs.}
  
  \item{corr}{Data frame of correlation results between NES values of 
  different inputs.}
  
  \item{corr.plot}{Tile plot of correlation results between NES values of 
  different inputs.}
}

\author{Belinda B. Garana
}

\note{Output venn.diagram will only display results for the first 4 data frames 
in the input data.list.
}

\seealso{
  \link[panSEA]{mGSEA}
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
weights$Rank_metric <- rnorm(length(Gene))
weights2$Rank_metric <- rnorm(length(Gene))

## Step 2: prepare set info
info <- as.data.frame(Gene)

# moa is our default name for the column containing drug set annotations
info$gs_name <- rep(paste("Set", LETTERS[seq(from = 1, to = 5)]), 10)

# convert data frame into gmt object
gmt <- DMEA::as_gmt(info, element.names = "Gene", set.names = "gs_name")

## Step 3: perform ssGSEA and store results
ssGSEA.results <- list()
ssGSEA.results[["Transcriptomics"]] <- panSEA::ssGSEA(weights, gmt)
ssGSEA.results[["Proteomics"]] <- panSEA::ssGSEA(weights2, gmt)

## Step 4: compile list of ssGSEA results
compile_mGSEA.test <- panSEA::compile_mGSEA(ssGSEA.results)
}
