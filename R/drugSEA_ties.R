
# load GSEA_custom function
GSEA_custom <- function(input.df, gmt.list,
                        nperm = 1000,
                        stat.type = "Weighted", min.per.set,
                        convert.synonyms = FALSE) {
  GSEA.EnrichmentScore <- function(gene.list, gene.set,
                                   weighted.score.type = score.weight,
                                   correl.vector = NULL) {
    tag.indicators <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicators
    N <- length(gene.list)
    Nh <- sum(tag.indicators)
    Nm <- N - Nh
    if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicators == 1])
    norm.tag <- 1.0 / sum.correl.tag
    norm.no.tag <- 1.0 / Nm
    RES <-
      cumsum(tag.indicators * correl.vector * norm.tag -
               no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(
      ES = ES,
      arg.ES = arg.ES,
      RES = RES,
      indicator = tag.indicators
    ))
  }
  
  # set weight for drug ranks
  if (stat.type == "Classic") {
    score.weight <- 0
  } else if (stat.type == "Weighted") {
    score.weight <- 1
  }
  
  # find biggest drug set
  biggest.set <- min.per.set
  for (i in seq_len(length(gmt.list$genesets))) {
    if (length(gmt.list$genesets[[i]]) > biggest.set) {
      biggest.set <- length(gmt.list$genesets[[i]])
    }
  }
  
  # restructure drug sets so that each MOA is a column with drugs in rows
  # also generate list of all annotated drugs from gmt
  Drug.Sets <- matrix(
    data = "", nrow = biggest.set,
    ncol = length(gmt.list$genesets)
  )
  colnames(Drug.Sets) <- gmt.list$geneset.names
  annotated.drugs <- c()
  for (i in seq_len(ncol(Drug.Sets))) {
    for (j in seq_len(length(gmt.list$genesets[[i]]))) {
      Drug.Sets[j, i] <- gmt.list$genesets[[i]][[j]][1]
      annotated.drugs <- c(annotated.drugs, gmt.list$genesets[[i]][[j]][1])
    }
  }
  Drug.Sets <- as.data.frame(Drug.Sets)
  annotated.drugs <- unique(annotated.drugs)
  
  # check class of inputs is as expected
  testthat::expect_is(input.df, "data.frame")
  testthat::expect_is(Drug.Sets, "data.frame")
  
  Samples <- colnames(input.df)[2]
  Drug.Sets.All <- colnames(Drug.Sets)
  
  annotations <- matrix(
    data = 0, nrow = nrow(input.df),
    ncol = length(Drug.Sets.All)
  )
  colnames(annotations) <- Drug.Sets.All
  
  annotations <- as.data.frame(annotations)
  
  annotations <- cbind(input.df[, 1], annotations)
  colnames(annotations) <- c(colnames(input.df)[1], Drug.Sets.All)
  annotations <- as.matrix(annotations)
  
  # Find out how many cores are available
  if (requireNamespace("parallel") &
      requireNamespace("snow") &
      requireNamespace("doSNOW")) {
    cores <- parallel::detectCores()
    if (cores[1] > 1) {
      cl <- snow::makeCluster(cores[1] - 1) # cluster using all but 1 core
      doSNOW::registerDoSNOW(cl) # register cluster
    }
  }
  
  num.hits.pathways <- as.data.frame(Drug.Sets.All)
  num.hits.pathways$N_drugs <- NA
  ### Annotate drug sets
  for (j in seq_len(length(Drug.Sets.All))) {
    temp.pathway <- Drug.Sets[, Drug.Sets.All[j]]
    for (i in seq_len(nrow(annotations))) {
      if (annotations[i, 1] %in% temp.pathway) {
        annotations[i, j + 1] <- "X"
      }
    }
    num.hits.pathways$N_drugs[j] <-
      sum(annotations[, Drug.Sets.All[j]] == "X")
  }
  
  ### Check if any drugs are not in any sets
  unannotated.drugs <- c()
  for (i in seq_len(nrow(annotations))) {
    if (sum(annotations[i, 2:ncol(annotations)] == "X") == 0) {
      unannotated.drugs <- c(unannotated.drugs, annotations[i, 1])
    }
  }
  
  ### Check if these unannotated drug names are a synonym for an annotated drug
  if (convert.synonyms) {
    # access drug synonyms
    synonyms <- utils::read.csv(
      paste0("https://raw.github.com/BelindaBGarana/DMEA/",
             "shiny-app/Inputs/PRISM_drug_synonyms.csv")
    )
    replacements <- as.data.frame(unannotated.drugs)
    replacements[, c("PRISM_drug_name",
                     "PubChem_CID", "PubChem_synonyms")] <- NA
    for (i in seq_len(length(unannotated.drugs))) {
      # check if unannotated drug name is a CID for a PRISM drug
      if (unannotated.drugs[i] %in% synonyms$PubChem_CID) {
        replacements[i, c("PRISM_drug_name",
                          "PubChem_CID", "PubChem_synonyms")] <-
          synonyms[synonyms$PubChem_CID == unannotated.drugs[i], ]
      } else {
        # else check if unannotated drug is a synonym for a PRISM drug name
        for (j in seq_len(nrow(synonyms))) {
          temp.synonyms <- strsplit(synonyms$PubChem_synonyms[j], "[|]")
          if (unannotated.drugs[i] %in% temp.synonyms) {
            replacements[i, c("PRISM_drug_name",
                              "PubChem_CID", "PubChem_synonyms")] <-
              synonyms[j, c("PRISM_drug_name",
                            "PubChem_CID", "PubChem_synonyms")]
            break
          }
        }
      }
    }
    unannotated.drugs <-
      replacements[is.na(replacements$PRISM_drug_name), ]$unannotated.drugs
    
    # remove any drugs which were not replaced
    replacements <- na.omit(replacements)
    rm(synonyms)
    
    # replace drug names in rank list
    for (i in seq_len(nrow(input.df))) {
      if (input.df[i, 1] %in% replacements$unannotated.drugs) {
        input.df[i, 1] <-
          replacements[replacements$unannotated.drugs == input.df[i, 1],
          ]$PRISM_drug_name
      }
    }
    
    # redo annotations
    annotations <- matrix(
      data = 0, nrow = nrow(input.df),
      ncol = length(Drug.Sets.All)
    )
    colnames(annotations) <- Drug.Sets.All
    
    annotations <- as.data.frame(annotations)
    
    annotations <- cbind(input.df[, 1], annotations)
    colnames(annotations) <- c(colnames(input.df)[1], Drug.Sets.All)
    annotations <- as.matrix(annotations)
    
    ### Annotate drug sets
    for (j in seq_len(length(Drug.Sets.All))) {
      temp.pathway <- Drug.Sets[, Drug.Sets.All[j]]
      for (i in seq_len(nrow(annotations))) {
        if (annotations[i, 1] %in% temp.pathway) {
          annotations[i, j + 1] <- "X"
        }
      }
      num.hits.pathways$N_drugs[j] <-
        sum(annotations[, Drug.Sets.All[j]] == "X")
    }
  } else {
    replacements <- NULL
  }
  unannotated.input <- input.df[input.df[, 1] %in% unannotated.drugs, ]
  
  Drug.Sets.under <-
    num.hits.pathways[num.hits.pathways$N_drugs < min.per.set, ]$Drug.Sets.All
  Drug.Sets.over <-
    num.hits.pathways[num.hits.pathways$N_drugs >= min.per.set,
    ]$Drug.Sets.All
  if (length(Drug.Sets.under) > 0) {
    warning(paste0(
      "Removing drug sets with less than ", min.per.set,
      " drugs observed in data set"
    ))
    annotations <- annotations[, c(colnames(input.df)[1], Drug.Sets.All)]
    
    # only keep drug sets with >= min.per.set drugs
    Drug.Sets.All <- Drug.Sets.over
  }
  colnames(num.hits.pathways)[1] <- "Drug_set"
  if (length(Drug.Sets.over) < 2) {
    stop("annotations for 2+ drug sets are required")
  } # BG 2022-06-25
  annotations <- as.data.frame(annotations)
  input.df <- merge(input.df, annotations, by = colnames(input.df)[1])
  input.df <- stats::na.omit(input.df)
  rm(annotations)
  
  GSEA.Results <- matrix(data = NA, nrow = length(Drug.Sets.All), ncol = 8)
  colnames(GSEA.Results) <- c(
    "Rank_metric", "Drug_set", "ES", "NES",
    "p_value", "Position_at_max", "FDR_q_value", "Leading_edge"
  )
  GSEA.Results <- as.data.frame(GSEA.Results)
  GSEA.Results$Drug_set <- Drug.Sets.All
  GSEA.Results$Rank_metric <- Samples
  
  ## Assuming first two columns in data table are drug names and
  ## rank metric (e.g. Foldchange, SNR)
  input.df[, Samples] <- as.numeric(as.character(input.df[, Samples]))
  
  # sort by descending order for the rank metric
  input.df <- input.df[order(-input.df[, Samples]), ]
  
  # reorder row indices for counting in for loop below
  rownames(input.df) <- seq_len(nrow(input.df))
  ions <- nrow(input.df)
  
  # for plotting
  ks_results_plot <- list()
  positions.of.hits <- list()
  
  gene.list <- seq_len(ions)
  rank_metric <- input.df[, Samples] # Save the rank metric
  
  pos_gene_set <- array(data = 0, dim = nrow(input.df), dimnames = NULL)
  
  ## Calculate Real KS Statistic with original input
  message("calculating ES with original input")
  for (i in seq_len(length(Drug.Sets.All))) {
    pos_gene_set <- which(input.df[, Drug.Sets.All[i]] %in% c("X"))
    KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set,
                                    weighted.score.type = score.weight,
                                    correl.vector = rank_metric
    )
    
    # store ES and positions
    GSEA.Results[GSEA.Results$Drug_set == Drug.Sets.All[i],
    ]$ES <- KS_real$ES
    GSEA.Results[GSEA.Results$Drug_set == Drug.Sets.All[i],
    ]$Position_at_max <- KS_real$arg.ES
    ks_results_plot[[Drug.Sets.All[i]]] <- KS_real$RES
    positions.of.hits[[Drug.Sets.All[i]]] <- pos_gene_set
    
    # identify leading edge elements
    if (KS_real$ES <= 0) {
      pos_leading_edge <- pos_gene_set[pos_gene_set >= KS_real$arg.ES]
    } else if (KS_real$ES >= 0) {
      pos_leading_edge <- pos_gene_set[pos_gene_set <= KS_real$arg.ES]
    }
    leading_edge <- input.df[pos_leading_edge, 1]
    GSEA.Results[GSEA.Results$Drug_set == Drug.Sets.All[i], 
    ]$Leading_edge <- paste(leading_edge, collapse = ", ")
  }
  Mountain.Plot.Info <- list(
    MountainPlot = ks_results_plot,
    Position.of.hits = positions.of.hits
  )
  rm(pos_gene_set)
  rm(KS_real)
  
  ### tie permutations (NOT a null distribution)
  KSRandomArray_ties <- matrix(data = NA, nrow = nperm,
                               ncol = length(Drug.Sets.All))
  num.Drug.Sets.all <- length(Drug.Sets.All)
  `%dopar%` <- foreach::`%dopar%`
  message("calculating ES with shuffled ties")
  KSRandomArray_ties <- foreach::foreach(L = seq_len(nperm),
                                         .combine = "rbind") %dopar% {
                                           # randomly jitter list to shuffle ties
                                           temp.input <- input.df[order(-input.df[, Samples], jitter(-input.df[, Samples])),]
                                           rownames(temp.input) <- seq_len(nrow(temp.input)) # adding in case random result is due to nonsequential rownames?
                                           gene.list <- rownames(temp.input)
                                           
                                           temp.KSRandomArray <- matrix(data = NA, nrow = 1,
                                                                        ncol = num.Drug.Sets.all)
                                           for (i in seq_len(length(Drug.Sets.All))) {
                                             pos_gene_set <- which(temp.input[, Drug.Sets.All[i]] %in% c("X"))
                                             
                                             # calculate ES with shuffled ties
                                             temp.KSRandomArray[, i] <- GSEA.EnrichmentScore(gene.list,
                                                                                             pos_gene_set,
                                                                                             weighted.score.type = score.weight,
                                                                                             correl.vector = temp.input[, Samples]
                                             )$ES
                                           }
                                           temp.KSRandomArray
                                         }
  colnames(KSRandomArray_ties) <- Drug.Sets.All
  KSRandomArray_ties <- as.data.frame(KSRandomArray_ties)
  
  ### calculate ES & NES
  # let ES = mean of ES tie distribution and NES = mean of NES tie distribution
  GSEA.Results.ties <- matrix(data = NA, nrow = length(Drug.Sets.All), ncol = 11)
  colnames(GSEA.Results.ties) <- c(
    "Rank_metric", "Drug_set", "ES","ES_min","ES_max", "ES_sd",
    "NES", "NES_min", "NES_max", "p_value", "FDR_q_value"
  )
  GSEA.Results.ties <- as.data.frame(GSEA.Results.ties)
  GSEA.Results.ties$Drug_set <- Drug.Sets.All
  GSEA.Results.ties$Rank_metric <- Samples
  GSEA.Results.ties$ES <- colMeans(KSRandomArray_ties, na.rm = TRUE)
  GSEA.Results.ties$ES_min <- apply(KSRandomArray_ties, 2, min)
  GSEA.Results.ties$ES_max <- apply(KSRandomArray_ties, 2, max)
  GSEA.Results.ties$ES_sd <- apply(KSRandomArray_ties, 2, sd)
  percent.pos.GSEA.ties <- sum(GSEA.Results.ties$ES >= 0) /
    length(GSEA.Results.ties$ES) # BG OR EQUAL TO
  percent.neg.GSEA.ties <- sum(GSEA.Results.ties$ES <= 0) /
    length(GSEA.Results.ties$ES) # BG OR EQUAL TO
  
  ### MOA permutations (null distribution)
  KSRandomArray <- matrix(data = NA, nrow = nperm,
                          ncol = length(Drug.Sets.All))
  num.Drug.Sets.all <- length(Drug.Sets.All)
  `%dopar%` <- foreach::`%dopar%`
  
  KSRandomArray <- foreach::foreach(L = seq_len(nperm),
                                    .combine = "rbind") %dopar% {
                                      temp.KSRandomArray <- matrix(data = NA, nrow = 1,
                                                                   ncol = num.Drug.Sets.all)
                                      for (i in seq_len(length(Drug.Sets.All))) {
                                        # randomly sample drug set of same size as original set
                                        numhits_pathway <- length(positions.of.hits[[Drug.Sets.All[i]]])
                                        pos_gene_set <- sample(seq_len(ions), numhits_pathway)
                                        
                                        # calculate null enrichment score
                                        temp.KSRandomArray[, i] <- GSEA.EnrichmentScore(gene.list,
                                                                                         pos_gene_set,
                                                                                         weighted.score.type = score.weight,
                                                                                         correl.vector = rank_metric
                                        )$ES
                                      }
                                      temp.KSRandomArray
                                    }
  colnames(KSRandomArray) <- Drug.Sets.All
  KSRandomArray <- as.data.frame(KSRandomArray)
  
  ### normalize the GSEA distribution
  KSRandomArray.Norm <- matrix(
    data = NA, nrow = nrow(KSRandomArray),
    ncol = ncol(KSRandomArray)
  )
  colnames(KSRandomArray.Norm) <- colnames(KSRandomArray)
  KSRandomArray.temp <- 0
  for (i in seq_len(ncol(KSRandomArray.Norm))) {
    KSRandomArray.temp <- KSRandomArray[, i]
    pos.temp <- KSRandomArray.temp[which(KSRandomArray.temp >= 0)]
    neg.temp <-
      KSRandomArray.temp[which(KSRandomArray.temp <= 0)] # BG OR EQUAL TO
    
    avg.pos <- mean(pos.temp)
    avg.neg <- mean(neg.temp)
    
    norm.pos.temp <- pos.temp / avg.pos
    norm.neg.temp <- neg.temp / avg.neg * -1
    
    norm.perms <- c(norm.pos.temp, norm.neg.temp)
    
    KSRandomArray.Norm[, i] <- norm.perms
  }
  GSEA.NES.perms <- as.vector(KSRandomArray.Norm)
  rm(KSRandomArray.Norm)
  GSEA.NES.perms.pos <- GSEA.NES.perms[which(GSEA.NES.perms >= 0)]
  GSEA.NES.perms.neg <-
    GSEA.NES.perms[which(GSEA.NES.perms <= 0)] # BG OR EQUAL TO
  rm(GSEA.NES.perms)
  percent.pos.GSEA <- sum(GSEA.Results$ES >= 0) /
    length(GSEA.Results$ES) # BG OR EQUAL TO
  percent.neg.GSEA <- sum(GSEA.Results$ES <= 0) /
    length(GSEA.Results$ES) # BG OR EQUAL TO
  
  # Calculate GSEA NES and p-value
  for (i in seq_len(length(Drug.Sets.All))) {
    temp.gene.set <- Drug.Sets.All[i]
    temp.perms <- KSRandomArray[, temp.gene.set]
    temp.ES <- GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$ES
    temp.ES.tie <- GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$ES
    temp.ES.tie.min <- GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$ES_min
    temp.ES.tie.max <- GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$ES_max
    pos.perms <- temp.perms[which(temp.perms >= 0)] # BG OR EQUAL TO
    neg.perms <- temp.perms[which(temp.perms <= 0)] # BG OR EQUAL TO
    
    # original input order
    if (temp.ES >= 0) { # BG OR EQUAL TO
      # p-value
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$p_value <-
        signif(sum(pos.perms >= temp.ES) /
                 length(pos.perms), digits = 3) # BG OR EQUAL TO
      
      # NES
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$NES <-
        signif(temp.ES / mean(pos.perms), digits = 3)
    } else { # BG OR EQUAL TO
      # p-value
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$p_value <-
        signif(sum(neg.perms <= temp.ES) /
                 length(neg.perms), digits = 3) # BG OR EQUAL TO
      
      # NES
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$NES <-
        signif(temp.ES / mean(neg.perms) * -1, digits = 3)
    }
    
    # shuffled ties
    if (temp.ES.tie >= 0) {
      # p-value
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$p_value <-
        signif(sum(pos.perms >= temp.ES.tie.min) /
                 length(pos.perms), digits = 3)
      
      # NES
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES <-
        signif(temp.ES.tie / mean(pos.perms), digits = 3)
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES_min <-
        signif(temp.ES.tie.min / mean(pos.perms), digits = 3)
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES_max <-
        signif(temp.ES.tie.max / mean(pos.perms), digits = 3)
    } else { # BG OR EQUAL TO
      # p-val
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$p_value <-
        signif(sum(neg.perms <= temp.ES) /
                 length(neg.perms), digits = 3) # BG OR EQUAL TO
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$p_value <-
        signif(sum(neg.perms <= temp.ES.tie.max) /
                 length(neg.perms), digits = 3)
      
      # NES
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES <-
        signif(temp.ES.tie / mean(neg.perms) * -1, digits = 3)
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES_min <-
        signif(temp.ES.tie.min / mean(neg.perms) * -1, digits = 3)
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES_max <-
        signif(temp.ES.tie.max / mean(neg.perms) * -1, digits = 3)
    }
      
  }
  
  # Calculate GSEA FDR
  for (i in seq_len(length(Drug.Sets.All))) {
    temp.gene.set <- Drug.Sets.All[i]
    temp.NES <- GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$NES
    temp.NES.tie <- GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES
    temp.NES.tie.min <- GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES_min
    temp.NES.tie.max <- GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$NES_max
    
    # results for original input
    if (temp.NES >= 0) { # BG OR EQUAL TO
      percent.temp <-
        sum(GSEA.NES.perms.pos >=
              GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$NES) /
        length(GSEA.NES.perms.pos) # BG OR EQUAL TO
      percent.pos.stronger <-
        sum(GSEA.Results$NES >=
              GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$NES) /
        sum(GSEA.Results$NES >= 0) # BG
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$FDR_q_value <-
        ifelse(signif(percent.temp / percent.pos.stronger, digits = 3) < 1,
               signif(percent.temp / percent.pos.stronger, digits = 3), 1
        ) # BG
    } else { # BG OR EQUAL TO
      percent.temp <- sum(GSEA.NES.perms.neg <= temp.NES) /
        length(GSEA.NES.perms.neg) # BG OR EQUAL TO
      percent.neg.stronger <- sum(GSEA.Results$NES <= temp.NES) /
        sum(GSEA.Results$NES <= 0) # BG
      GSEA.Results[GSEA.Results$Drug_set == temp.gene.set, ]$FDR_q_value <-
        ifelse(signif(percent.temp / percent.neg.stronger, digits = 3) < 1,
               signif(percent.temp / percent.neg.stronger, digits = 3), 1
        ) # BG
    }
    
    # results for shuffled ties
    if (temp.NES.tie >= 0) {
      percent.temp <-
        sum(GSEA.NES.perms.pos >= temp.NES.tie.min) / length(GSEA.NES.perms.pos) # BG OR EQUAL TO
      percent.pos.stronger <-
        sum(GSEA.Results.ties$NES_max >= temp.NES.tie.min) / sum(GSEA.Results.ties$NES >= 0) # BG
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$FDR_q_value <-
        ifelse(signif(percent.temp / percent.pos.stronger, digits = 3) < 1,
               signif(percent.temp / percent.pos.stronger, digits = 3), 1
        ) # BG
    } else {
      percent.temp <- sum(GSEA.NES.perms.neg <= temp.NES.tie.max) / length(GSEA.NES.perms.neg) # BG OR EQUAL TO
      percent.neg.stronger <- sum(GSEA.Results.ties$NES_min <= temp.NES.tie.max) / sum(GSEA.Results.ties$NES <= 0) # BG
      GSEA.Results.ties[GSEA.Results.ties$Drug_set == temp.gene.set, ]$FDR_q_value <-
        ifelse(signif(percent.temp / percent.neg.stronger, digits = 3) < 1,
               signif(percent.temp / percent.neg.stronger, digits = 3), 1
        ) # BG
    }
  }
  
  if (requireNamespace("parallel") &
      requireNamespace("snow") &
      requireNamespace("doSNOW")) {
    if (cores[1] > 1) {
      snow::stopCluster(cl) # stop cluster
      rm(cl)
    }
  }
  
  return(list(
    GSEA.Results = GSEA.Results,
    GSEA.Results.ties = GSEA.Results.ties,
    Mountain.Plot.Info = Mountain.Plot.Info,
    ranking.metric = rank_metric,
    removed.sets =
      num.hits.pathways[num.hits.pathways$N_drugs < min.per.set, ],
    replacements = replacements,
    unannotated.drugs = unannotated.input
  ))
}

# load mountain plot function
gsea_mountain_plot <- function(GSEA.list, Sample.Name, Gene.Set.A,
                               color = TRUE) {
  if (color) {
    color.palette <- c("red")
  } else if (!color) {
    color.palette <- c("black")
  }
  
  GSEA.list <- GSEA.list
  plotting.info <- GSEA.list$Mountain.Plot.Info
  ES_profiles <- plotting.info$MountainPlot
  position.of.hits <- plotting.info$Position.of.hits
  
  ranking.metric <- as.numeric(GSEA.list$ranking.metric)
  
  ES_GeneSetA <- ES_profiles[[Gene.Set.A]]
  
  pos.GeneSetA <- position.of.hits[[Gene.Set.A]]
  
  myName <- Gene.Set.A
  
  GSEA.results <- GSEA.list$GSEA.Results
  GSEA.results <- GSEA.results[GSEA.results$Rank_metric == Sample.Name, ]
  GSEA.results$p_value <- as.numeric(as.character(GSEA.results$p_value))
  GSEA.results$FDR_q_value <- as.numeric(as.character(GSEA.results$FDR_q_value))
  
  gsea.normalized.enrichment.score <-
    signif(GSEA.results[GSEA.results$Drug_set == myName, ]$NES, digits = 3)
  gsea.p.value <-
    signif(GSEA.results[GSEA.results$Drug_set == myName, ]$p_value,
           digits = 3)
  gsea.q.val <-
    signif(GSEA.results[GSEA.results$Drug_set == myName, ]$FDR_q_value,
           digits = 3)
  
  results.as.text <- paste(
    "NES", gsea.normalized.enrichment.score, "\np-value:", gsea.p.value,
    "\nq-value:", gsea.q.val
  )
  
  gsea.hit.indices.A <- pos.GeneSetA
  
  gsea.es.profile.A <- ES_GeneSetA[gsea.hit.indices.A]
  
  ngenes <- length(ES_GeneSetA)
  combined.ES <- c(ES_GeneSetA)
  max.ES <- max(combined.ES)
  min.ES <- min(combined.ES)
  
  mtn.plot <- NULL
  mtn.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(
        x = c(0, pos.GeneSetA, ngenes),
        y = c(0, gsea.es.profile.A, 0)
      ),
      color = color.palette[1], size = 1
    ) +
    ggplot2::labs(x = NULL, y = "ES", title = myName) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(0, "pt"),
      axis.title.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  if (!is.null(gsea.normalized.enrichment.score) &
      length(gsea.normalized.enrichment.score) > 0) {
    if (gsea.normalized.enrichment.score > 0) {
      mtn.plot <- mtn.plot + ggplot2::annotate("text",
                                               label = results.as.text, x = ngenes * 0.90,
                                               y = max.ES * 0.70, color = "black"
      )
    } else if (gsea.normalized.enrichment.score < 0) {
      mtn.plot <- mtn.plot + ggplot2::annotate("text",
                                               label = results.as.text, x = ngenes * 0.10,
                                               y = min.ES * 0.30, color = "black"
      )
    }
  }
  
  hit.markers <- ggplot2::ggplot() +
    cowplot::theme_nothing() +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black",
                                                        fill = NA)) +
    ggplot2::geom_col(ggplot2::aes(x = pos.GeneSetA, y = 1),
                      color = color.palette[1],
                      fill = color.palette[1], width = 0.75
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    aplot::xlim2(mtn.plot)
  
  show.rank.metric <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(ngenes),
      y = ranking.metric
    ), color = "grey50") +
    ggplot2::geom_area(ggplot2::aes(
      x = seq_len(ngenes),
      y = ranking.metric
    ), fill = "grey50") +
    ggplot2::labs(x = "Rank", y = "Rank Metric") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
    ) +
    aplot::xlim2(mtn.plot)
  
  assembled <- mtn.plot + hit.markers + show.rank.metric +
    patchwork::plot_layout(ncol = 1, heights = c(2.5, 0.5, 2))
  return(assembled)
}

summaryPlots <- function(EA, FDR = 0.25, n.top = 10, est.name = "Pearson.est") {
  plot.data <- EA$GSEA.Results.ties
  
  if (nrow(plot.data[plot.data$p_value == 0, ]) > 0) {
    plot.data[plot.data$p_value == 0, ]$p_value <- 0.00099
  }
  
  # define x, y limits
  limit.x <- ceiling(max(abs(as.numeric(plot.data$NES)), na.rm = TRUE))
  limit.y <- ceiling(max(-as.numeric(log(plot.data$p_value, 10)), na.rm = TRUE))
  
  ## select results which meet FDR threshold and produce mountain plots
  significant.hits <- EA$GSEA.Results.ties[which(EA$GSEA.Results.ties$FDR_q_value < FDR), ]
  temp.plot <- list()
  if (nrow(significant.hits) > 0) {
    for (i in seq_len(nrow(significant.hits))) {
      temp <- gsea_mountain_plot(
        GSEA.list = EA, Sample.Name = est.name,
        Gene.Set.A = significant.hits$Drug_set[i]
      )
      temp.plot[significant.hits$Drug_set[i]] <- list(temp)
    }
  } else {
    warning("No enrichments met the FDR cut-off to produce mountain plots")
  }
  
  # categorize data by significance level if there are significant hits
  if (nrow(significant.hits) > 0) {
    plot.data$Significance <- paste0("FDR > ", FDR)
    plot.data[plot.data$FDR_q_value < FDR, ]$Significance <- paste0("FDR < ", FDR)
    plot.data$Significance <- factor(plot.data$Significance,
                                     levels = c(paste0("FDR < ", FDR),
                                                paste0("FDR > ", FDR)))
    bar.data <- plot.data %>% slice_max(abs(NES), n=n.top)
    bar.data <- bar.data[order(bar.data$NES),]
    
    volc <- ggplot2::ggplot(data = plot.data, aes(
      x = NES, y = -log(p_value, 10),
      color = Significance
    )) +
      ggplot2::geom_point(size = 4) +
      ggrepel::geom_text_repel(
        data = subset(plot.data, Significance == paste0("FDR < ", FDR)),
        mapping = aes(label = Drug_set, size = I(6)), nudge_y = 0.25
      ) +
      ggplot2::scale_color_manual(
        values = c("red", "azure4"), name = "Significance",
        breaks = c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))
      ) +
      ggplot2::xlim(-limit.x, limit.x) +
      ggplot2::ylim(0, limit.y) +
      ggplot2::xlab("Normalized Enrichment Score") +
      ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                          linewidth = 0.5) +
      ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 26, face = "bold"),
        panel.background = element_rect(
          fill = "white", colour = "white", linewidth = 0.5,
          linetype = "solid", color = "black"
        ), text = element_text(size = 20),
        legend.position = "bottom", legend.key = element_blank()
      )
    
    bar <- ggplot2::ggplot(bar.data, aes(x=NES, y=Drug_set, fill = Significance)) + 
      geom_bar(stat='identity') + ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 26, face = "bold"),
        panel.background = element_rect(
          fill = "white", colour = "white", linewidth = 0.5,
          linetype = "solid", color = "black"
        ), text = element_text(size = 20),
        legend.position = "bottom", legend.key = element_blank()
      ) + ggplot2::xlab("Normalized Enrichment Score") + 
      ggplot2::ylab("") +
      ggplot2::scale_fill_manual(
        values = c("red", "azure4"), name = "Significance",
        breaks = c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))
      ) 
  } else {
    plot.data$Significance <- paste0("FDR > ", FDR)
    plot.data$Significance <- factor(plot.data$Significance,
                                     levels = paste0("FDR > ", FDR))
    bar.data <- plot.data %>% slice_max(abs(NES), n=n.top)
    bar.data <- bar.data[order(bar.data$NES),]
    
    volc <- ggplot2::ggplot(data = plot.data,
                            aes(x = NES, y = -log(p_value, 10))) +
      ggplot2::geom_point(size = 4, color = "azure4") +
      ggplot2::xlim(-limit.x, limit.x) +
      ggplot2::ylim(0, limit.y) +
      ggplot2::xlab("Normalized Enrichment Score") +
      ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                          linewidth = 0.5) +
      ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 26, face = "bold"),
        panel.background = element_rect(
          fill = "white", colour = "white", linewidth = 0.5,
          linetype = "solid", color = "black"
        ), text = element_text(size = 20),
        legend.position = "bottom", legend.key = element_blank()
      )
    
    bar <- ggplot2::ggplot(bar.data, aes(x=NES, y=Drug_set, fill = Significance)) + 
      geom_bar(stat='identity') + ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 26, face = "bold"),
        panel.background = element_rect(
          fill = "white", colour = "white", linewidth = 0.5,
          linetype = "solid", color = "black"
        ), text = element_text(size = 20),
        legend.position = "bottom", legend.key = element_blank()
      ) + ggplot2::xlab("Normalized Enrichment Score") + 
      ggplot2::ylab("") +
      ggplot2::scale_fill_manual(
        values = c("azure4"), name = "Significance",
        breaks = c(paste0("FDR > ", FDR))
      ) 
  }
  bar.data$name <- est.name
  bg.theme <- ggplot2::theme(
    legend.background = element_rect(), legend.position = "top",
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    legend.title = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(
      lineheight = .8, face = "bold", size = 36
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )
  dot <- ggplot2::ggplot(
    bar.data,
    ggplot2::aes(x = name, y = Drug_set, color = NES,
      size = -log10(FDR_q_value)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = bar.data[order(bar.data$NES, decreasing = TRUE), ]$Drug_set) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      y = "",
      color = "NES", size = "-log(FDR)"
    )
  
  dot.sd <- ggplot2::ggplot(
    bar.data,
    ggplot2::aes(x = name, y = Drug_set, color = NES,
                 size = -log10(ES_sd)
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = bar.data[order(bar.data$NES, decreasing = TRUE), ]$Drug_set) +
    viridis::scale_color_viridis() +
    bg.theme +
    ggplot2::labs(
      y = "",
      color = "NES", size = "-log(Standard Deviation of ES)"
    )
  return(list(volcano = volc, bar = bar, mtn = temp.plot, dot = dot, dot.sd = dot.sd))
}

drugSEA_ties <- function(data, gmt = NULL, drug = "Drug",
                    rank.metric = "Pearson.est", set.type = "moa",
                    direction.adjust = NULL, FDR = 0.25,
                    num.permutations = 1000, stat.type = "Weighted",
                    min.per.set = 6, n.top = 10, sep = "[|]",
                    exclusions = c("-666", "NA", "na", "NaN", "NULL"),
                    descriptions = NULL, convert.synonyms = FALSE, plots = TRUE) {
  # avoid unclear initialization
  Drug_set <- NULL
  NES <- NULL
  Significance <- NULL
  p_value <- NULL

  # check class of input is as expected
  testthat::expect_is(data, "data.frame")

  ## remove any drugs with rank.metric = 0, NA
  ## because signed ranks are required for enrichment analysis
  data[, c(rank.metric)] <- as.numeric(data[, c(rank.metric)])
  data <- data[data[, c(rank.metric)] != 0 & !is.na(data[, c(rank.metric)]), ]

  ## check for duplicate entries for each drug
  ## because only 1 rank metric is allowed per drug
  drug.names <- unique(data[, c(drug)])
  if (length(drug.names) < nrow(data)) {
    stop("each drug must have only 1 rank.metric value")
  }

  ## if necessary, generate gmt object
  if (is.null(gmt)) {
    gmt <- as_gmt(data,
      element.names = drug, set.names = set.type,
      min.per.set, sep, exclusions, descriptions
    )
  }

  ## if necessary, adjust rank metrics
  if (length(direction.adjust) > 0) {
    info <- dplyr::distinct(data[, c(drug, rank.metric, set.type)])
    info$adjust <- 1
    info[grep(direction.adjust, info[, c(set.type)],
              value = TRUE), "adjust"] <- -1
    info$adjusted.est <-
      as.numeric(info[, c(rank.metric)]) * as.numeric(info$adjust)
    input <- info[, c(drug, "adjusted.est")]
    est.name <- paste0(rank.metric, " direction-adjusted")
  } else {
    input <- data[, c(drug, rank.metric)]
    est.name <- rank.metric
  }
  colnames(input) <- c(drug, est.name)

  ## perform enrichment analysis using GSEA_custom
  message("Running enrichment analysis...")
  EA <- GSEA_custom(
    input, gmt, num.permutations, stat.type,
    min.per.set, convert.synonyms
  )
  
  if (plots) {
    ## produce volcano and bar plots
    sumPlots <- summaryPlots(EA, FDR = FDR, n.top = n.top, est.name = est.name) 
  } else {
    sumPlots <- list("volcano" = NA, "bar" = NA, "mtn" = NA)
  }

  return(list(
    gmt = gmt, result = EA$GSEA.Results.ties, result.w.ties = EA$GSEA.Results, 
    mtn.plots = sumPlots$mtn,
    volcano.plot = sumPlots$volcano, bar.plot = sumPlots$bar,
    dot.plot = sumPlots$dot, dot.sd = sumPlots$dot.sd,
    removed.sets = EA$removed.sets,
    replaced.drugs = EA$replacements,
    unannotated.drugs = EA$unannotated.drugs
  ))
}
