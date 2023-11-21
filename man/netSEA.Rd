\name{netSEA}

\alias{netSEA}

\title{Network Graph
}

\description{Generate network graph of leading edge elements 
from top significant enrichments
}

\usage{
netSEA(inputs, outputs, element.name = rep("Gene", length(inputs)), 
rank.var = rep("Log2FC", length(inputs)), 
n.top.sets = 3, p = 0.05, FDR = 0.25)
}

\arguments{
  \item{inputs}{List of data frames containing inputs used for 
  enrichment analyses. Each data frame must have column names matching the 
  element.name and rank.var given in the same order (e.g., if "Gene" and 
  "Log2FC" are the first element.name and rank.var, respectively, the first 
  data frame in the inputs list must contain columns with matching names).
  }

  \item{outputs}{List of data frames containing outputs from enrichment 
  analyses. Each data frame in the list must have columns "NES", "p-value", 
  "FDR_q_value", and "Leading_edge".
  }

  \item{element.name}{Vector of column names containing element names (e.g., 
  gene symbols) in each corresponding data frame in the inputs list.
  }

  \item{rank.var}{Vector of column names containing rank variable (e.g., 
  Log2FC) in each corresponding data frame in the inputs list.
  }

  \item{n.top.sets}{Number of sets from enrichment analyses to be represented 
  in the network graph. The top significant sets will be selected based on 
  abs(NES).
  }
  
  \item{p}{P-value cut-off for significant enrichments. Default: 0.05.
  }

  \item{FDR}{False discovery rate cut-off value for significant enrichments. 
  Default: 0.25.
  }
}

\value{Network graph of leading-edge elements from top significant enrichments.}

\author{Belinda B. Garana
}

\seealso{
  \link[netSEA]{compilDMEA}
}

\examples{
### Step 1: Run mGSEA (or mDMEA)
## Step 1a: prepare gene weight data frames
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

## Step 1b: prepare set info data frame
# alternatively, a gmt object could be provided
info <- as.data.frame(Gene)

# moa is our default name for the column containing drug set annotations
info$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 4)]), 6)

## Step 1c: perform DMEA and store results
data <- list(weights, weights2)
gmt <- DMEA::as_gmt(info)
types <- c("Transcriptomics", "Proteomics")
mGSEA.test <- mGSEA(data, gmt, types)

### Step 2: Extract results
outputs <- list()
for (i in 1:length(types)) {
  outputs[[types[i]]] <- ssGSEA.results[[types[i]]]$result
}

### Step 3: Create network graph
mGSEA.network <- netSEA(data, outputs)
}
