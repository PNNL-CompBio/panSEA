\name{mGSEA}

\alias{mGSEA}

\title{Multi-omic Gene Set Enrichment Analysis
}

\description{Runs enrichment analysis for each input omics type.
}

\usage{
mGSEA(data, gmt, types, direction.adjust=NULL, FDR=0.25, 
num.permutations=1000, stat.type="Weighted", min.per.set=6)
}

\arguments{
  \item{data}{List of data frames where the first column of each data frame 
  contains feature names (e.g., gene symbols) and the second column contains 
  the rank metric (e.g., Log2FC).
  }

  \item{gmt}{Gmt object containing sets of features (e.g., output from
  GSA.read.gmt function from GSA package).
  }

  \item{types}{Vector of character descriptions of each omics type (e.g., 
  c("Transcriptomics", "Proteomics")).
  }
  
  \item{direction.adjust}{Vector of conditions for set.type to identify which,
  if any, estimates to multiply by -1 (e.g., c(" agonist", " activator")).
  Default: NULL.
  }
  
  \item{FDR}{False discovery rate cut-off value to produce mountain plots for
  significant enrichments. Default: 0.25.
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
}

\value{
  \item{compiled.results}{Outputs from compilGSEA function.}

  \item{all.results}{List of outputs from ssGSEA function for each input type.}
}

\author{Belinda B. Garana
}

\section{Warning}{Each feature must have only one numeric, nonzero rank.
}

\seealso{
  \link[mGSEA]{compilDMEA}
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
info$gs_name <- rep(paste("Set", LETTERS[seq(from = 1, to = 4)]), 6)

# convert data frame into gmt object
gmt <- DMEA::as_gmt(info, element.names = "Gene", set.names = "gs_name")

## Step 3: perform DMEA and store results
mGSEA.test <- mGSEA(list(weights, weights2), gmt, 
                  types = c("Transcriptomics", "Proteomics"))
}
