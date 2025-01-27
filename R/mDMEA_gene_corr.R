mDMEA_gene_corr <- function(gmt = "PRISM",
                  expression = as.list(rep("L1000", length(types))),
                  weights, types, value = "Perturbation",
                  sample.names = colnames(expression[[1]])[1],
                  feature.names = rep("Gene", length(types)),
                  rank.metric = rep("Pearson.est", length(types)), p = 0.05,
                  FDR = 0.25, num.permutations = 1000, stat.type = "Weighted",
                  drug.info = NULL, drug = "Drug", set.type = "moa", 
                  min.per.set = 6, sep = "[|]", 
                  exclusions = c("-666", "NA", "na", "NaN", "NULL"),
                  descriptions = NULL, min.per.corr = 3, scatter.plots = TRUE,
                  scatter.plot.type = "pearson", FDR.scatter.plots = 0.05,
                  xlab = "Expression", ylab = value, 
                  position.x = "min", position.y = "min", se = TRUE, 
                  n.dot.sets = 10, ties = FALSEß) {
  #### Step 1. Load data if necessary ####
  # get gmt and expression if PRISM/L1000
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

  if ("L1000" %in% expression) {
    message("Loading L1000")
    # source: https://github.com/PNNL-CompBio/coderdata/blob/main/build/lincs/05-LINCS_perturbations.R
    # identify URLs
    #basename="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101406"
    L1000 <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5FL1000%5FLevel4%5FZSPCINF%5Fmlr12k%5Fn1667x12328.gctx.gz'
    L1000.genes <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5FL1000%5Fgene%5Finfo.txt.gz'
    L1000.inst <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101406/suppl/GSE101406%5FBroad%5FLINCS%5FL1000%5Finst%5Finfo.txt.gz'
    
    # download source sample & perturbation info
    L1000.genes.info <- readr::read_delim(L1000.genes, "\t") # pr_gene_symbol, pr_gene_id
    L1000.inst.info <- readr::read_delim(L1000.inst, "\t") # pert_iname; pert_type = "trt_cp" if drug, "ctl_vehicle" if DMSO; cell_id
    L1000.inst.info <- L1000.inst.info[L1000.inst.info$pert_type == "trt_cp", ] # only keep drug treatments (remove DMSO controls)
    
    # download L1000 data
    if (file.exists("L1000.gctx")) {
      L1000.df <- cmapR::parse_gctx("L1000.gctx")
    } else {
      res<-download.file(L1000,'L1000.gctx.gz', mode="wb")
      L1000.df<-cmapR::parse_gctx(R.utils::gunzip("L1000.gctx.gz")) 
    }
    
    # put data into long format
    L1000.long <- cmapR::melt_gct(L1000.df)
    colnames(L1000.long) <- c("inst_id", "pr_gene_id", "data_value")
    L1000.long$pr_gene_id <- as.numeric(L1000.long$pr_gene_id)
    
    # join with L1000.gene.info based on pr_gene_id;
    # L1000.inst.info based on inst_id;
    # keep only cell_id, pert_type, pert_iname, gene_symbol, data_value
    L1000.full <- L1000.long |>
      dplyr::left_join(L1000.genes.info) |>
      dplyr::left_join(L1000.inst.info) |>
      dplyr::select(cell_id,pert_type,pert_iname,pr_gene_symbol,data_value)|>
      dplyr::distinct()
    L1000.long <- NULL # save space
    
    # add relevant columns
    L1000.full <- na.omit(L1000.full)
    colnames(L1000.full)[1] <- "other_names" # match samples column name
    colnames(L1000.full)[4] <- "gene_symbol" # match genes column name
    colnames(L1000.full)[3] <- "chem_name"
    L1000.full$data_type <- "transcriptomics"
    L1000.full$source <- "Broad"
    L1000.full$study <- "LINCS"
    L1000.full$perturbation_type <- "drug" # all entries have pert_type="trt_cp"

    # aggregate across cell types
    L1000.full <- reshape2::dcast(L1000.full, gene_symbol ~ chem_name, 
                                  value.var = "data_value", fun.aggregate = mean)

    for (i in which(expression == "L1000")) {
      expression[[i]] <- L1000.full
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

    DMEA.list[[types[i]]] <- panSEA::DMEA_gene_corr(
      gmt, expression[[i]], weights[[i]], value, sample.names,
      feature.names[i], colnames(weights[[i]])[2],
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
