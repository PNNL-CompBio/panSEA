\name{compile_mDMEA}

\alias{compile_mDMEA}

\title{Compile Multi-omic Drug Mechanism Enrichment Analysis Results
}

\description{Compiles enrichment analysis results across different inputs.
}

\usage{
compile_mDMEA(DMEA.list, p=0.05, FDR=0.25, n.dot.sets=10)
}

\arguments{
  \item{DMEA.list}{List of outputs from different DMEA runs.
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

\references{
Garana, B.B., Joly, J.H., Delfarah, A. et al. Drug mechanism enrichment analysis
improves prioritization of therapeutics for repurposing. BMC Bioinformatics 24,
215 (2023). https://doi.org/10.1186/s12859-023-05343-8
}

\author{Belinda B. Garana
}

\note{Output venn.diagram will only display results for the first 4 data frames 
in the input data.list.
}

\seealso{
  \link[panSEA]{mDMEA}
}

\examples{
## Step 1: prepare drug weight data frames
# create list of drug symbols
Drug <- paste0("Drug_", seq(from = 1, to = 50))

# by default, drug names are found in
# the first column of your input data frame
ranks <- as.data.frame(Drug)
ranks2 <- as.data.frame(Drug)

# give each drug a rank
# by default, drug rank values are found in
# the second column of your input data frame
ranks$Pearson.est <- rnorm(length(Drug))
ranks2$Pearson.est <- rnorm(length(Drug))

## Step 2: prepare set info
info <- as.data.frame(Drug)

# moa is our default name for the column containing drug set annotations
info$gs_name <- rep(paste("Set", LETTERS[seq(from = 1, to = 5)]), 10)

# convert data frame into gmt object
gmt <- DMEA::as_gmt(info, element.names = "Drug", set.names = "gs_name")

## Step 3: perform drug set enrichment analysis and store results
DMEA.results <- list()
DMEA.results[["Transcriptomics"]] <- DMEA::drugSEA(ranks, gmt)
DMEA.results[["Proteomics"]] <- DMEA::drugSEA(ranks2, gmt)

## Step 4: compile list of DMEA results
compile_mDMEA.test <- panSEA::compile_mDMEA(DMEA.results)
}
