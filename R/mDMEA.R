mDMEA <- function(drug.sensitivity = "PRISM", gmt = "PRISM",
                  expression = as.list(rep("adherent CCLE", length(types))),
                  weights, types, value = "AUC",
                  sample.names = colnames(drug.sensitivity)[1],
                  feature.names = rep("Gene", length(types)),
                  weight.values = rep("Log2FC", length(types)),
                  rank.metric = rep("Pearson.est", length(types)), p = 0.05,
                  FDR = 0.25, num.permutations = 1000, stat.type = "Weighted",
                  drug.info = NULL, drug = "Drug", set.type = "moa", 
                  min.per.set = 6, sep = "[|]", 
                  exclusions = c("-666", "NA", "na", "NaN", "NULL"),
                  descriptions = NULL, min.per.corr = 3, scatter.plots = TRUE,
                  scatter.plot.type = "pearson", FDR.scatter.plots = 0.05,
                  xlab = "Weighted Voting Score", ylab = value, 
                  position.x = "min", position.y = "min", se = TRUE, 
                  n.dot.sets = 10, ties = FALSE) {
  #### Step 1. Load data if necessary ####
  # get drug.sensitivity, gmt, and expression if PRISM/CCLE
  if (is.character(drug.sensitivity)) {
    if (drug.sensitivity == "PRISM") {
      message("Loading PRISM drug sensitivity AUC scores")
      drug.sensitivity <- read.csv(file = paste0(
        "https://raw.github.com/BelindaBGarana/",
        "DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv"
      )) # 481 cell lines
      drug.sensitivity$X <- NULL
    } else {
      stop(paste(
        "drug.sensitivity must be either 'PRISM' or data frame per",
        "documentation"
      ))
    }
  }

  if (is.character(gmt)) {
    if (gmt == "PRISM") {
      message("Loading PRISM drug mechanism of action annotations")
      gmt <- GSA::GSA.read.gmt(
        file = paste0(
          "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
          "Inputs/MOA_gmt_file_n6_no_special_chars.gmt"
        )
      )
    } else {
      stop("gmt must be either 'PRISM' or list object per documentation")
    }
  }

  if ("adherent CCLE" %in% expression) {
    message("Loading adherent CCLE RNA-seq data version 19Q4")
    download.file(
      paste0(
        "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
      ),
      destfile =
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
    )
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")

    download.file(
      paste0(
        "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
      ),
      destfile =
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
    )
    load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")

    RNA.df <- rbind(RNA.first200, RNA.rest)

    for (i in which(expression == "adherent CCLE")) {
      expression[[i]] <- RNA.df
    }
  }

  # check that there are as many types as weights data frames
  if (length(types) != length(weights)) {
    stop("Length of types vector must match that of weights list")
  }

  # check that sample names are in drug.sensitivity data frame
  if (!(sample.names %in% names(drug.sensitivity))) {
    stop(paste("sample.names must match across drug.sensitivity and", 
               "expression data frames"))
  }

  # check that sample names are in expression data frames
  for (i in 1:length(types)) {
    if (!(sample.names %in% names(expression[[i]]))) {
      stop(paste("sample.names must match across drug.sensitivity and", 
                 "expression data frames"))
    }
  }

  #### Step 2. Perform DMEA on each omics type ####
  DMEA.list <- list()
  for (i in 1:length(types)) {
    message(paste("Running DMEA using", types[i], "data"))

    DMEA.list[[types[i]]] <- panSEA::DMEA(
      drug.sensitivity, gmt, expression[[i]],
      weights[[i]], value, sample.names,
      feature.names[i], weight.values[i],
      rank.metric[i], FDR, num.permutations,
      stat.type, drug.info, drug, set.type,
      min.per.set, sep, exclusions,
      descriptions, min.per.corr,
      scatter.plots, scatter.plot.type,
      FDR.scatter.plots, xlab, ylab,
      position.x, position.y, se, ties
    )

    # merge correlation results with drug annotations if !is.null(drug.info)
    if (!is.null(drug.info)) {
      DMEA.list[[types[i]]]$corr.result <-
        merge(DMEA.list[[types[i]]]$corr.result, drug.info,
          by.x = "Drug", by.y = drug
        )
    }
  }

  #### Step 3. Compile DMEA results across omics types ####
  if (length(types) > 1) {
    compiled.DMEA <- panSEA::compile_mDMEA(DMEA.list, p, FDR, n.dot.sets)
  } else {
    compiled.DMEA <- NA
  }

  return(list(
    compiled.results = compiled.DMEA,
    all.results = DMEA.list
  ))
}
