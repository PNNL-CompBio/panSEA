panSEA <- function(data, types, group.names="Diseased", 
                   group.samples=1:ncol(data[[1]]),
                   gmt.genes=NULL, gmt.drugs=NULL, gmt.ptm=NULL, 
                   ptm.split = "_", FDR=0.25, min.per.set=6, 
                   drug.sensitivity=NULL, expression=NULL){
  if (length(conditions) > 2 | length(conditions) < 1) {
    stop("Only 1 or 2 conditions are allowed")
  } else if (length(conditions) == 2) {
    deg <- list()
    deg.ptm <- list()
    #### Step 1. Differential expression analysis ####
    for (i in 1:length(types)) {
      if (grepl("phospho", types[i], ignore.case = TRUE)) {
        # separate data based on group (i.e., categorical phenotype)
        data1 <- data[[i]][ , group.samples[1]]
        data2 <- data[[i]][ , group.samples[2]]
        all.ptm.data <- cbind(data1, data2)
        
        # make sure there are no duplicated feature names
        all.ptm.data <- all.ptm.data[!duplicated(all.ptm.data), ]
        if (nrow(all.ptm.data) > unique(rownames(all.ptm.data))) {
          stop("There must be only one expression value per sample for each feature")
        }
        
        # create copy of data without phospho sites
        all.data <- all.ptm.data
        rownames(all.data) <- stringr::str_split_i(rownames(all.data),
                                                 ptm.split, 1)
        
        # with & without phospho-sites:
        # create expression sets
        eset.ptm <- Biobase::ExpressionSet(all.ptm.data)
        eset <- Biobase::ExpressionSet(all.data)
        
        # identify sample phenotypes
        gs <- factor(c(rep(0, ncol(data1)), rep(1, ncol(data2))))
        levels(gs) <- make.names(group.names)
        eset.ptm$group <- gs
        eset$group <- gs
        
        # set up design matrix
        design.ptm <- stats::model.matrix(~group + 0, eset.ptm)
        colnames(design.ptm) <- levels(gs)
        design <- stats::model.matrix(~group + 0, eset)
        colnames(design) <- levels(gs)
        
        # fit linear model
        fit.ptm <- limma::lmFit(eset.ptm, design.ptm)
        fit <- limma::lmFit(eset, design)
        
        # set up contrasts of interest and redo fit
        cts <- paste(groups[1], groups[2], sep = "-")
        cont.matrix.ptm <- limma::makeContrasts(contrasts=cts, levels=design.ptm)
        fit.ptm <- limma::contrasts.fit(fit.ptm, cont.matrix.ptm)
        cont.matrix <- limma::makeContrasts(contrasts=cts, levels=design)
        fit <- limma::contrasts.fit(fit, cont.matrix)
        
        # calculate differential expression and statistics
        fit.ptm <- limma::eBayes(fit.ptm, 0.01)
        deg.ptm[[types[i]]] <- limma::topTable(fit.ptm, adjust = "fdr", 
                                               number = nrow(fit))
        deg.ptm[[types[i]]]$Gene <- rownames(deg.ptm[[types[i]]])
        colnames(deg.ptm[[types[i]]])[1] <- "Log2FC"
        deg.ptm[[types[i]]] <- 
          deg.ptm[[types[i]]][ , c("Gene", 
                                   colnames(deg.ptm[[types[i]]])
                                   [1:(ncol(deg.ptm[[types[i]]])-1)])]
        
        fit <- limma::eBayes(fit, 0.01)
        deg[[types[i]]] <- limma::topTable(fit, adjust = "fdr", 
                                               number = nrow(fit))
        deg[[types[i]]]$Gene <- rownames(deg[[types[i]]])
        colnames(deg[[types[i]]])[1] <- "Log2FC"
        deg[[types[i]]] <- 
          deg[[types[i]]][ , c("Gene", 
                                   colnames(deg[[types[i]]])
                                   [1:(ncol(deg[[types[i]]])-1)])]


        
        tT <- limma::topTable(fit, adjust="fdr", sort.by="B", number=nrow(fit))
        tT$Gene <- rownames(tT)
        tT <- na.omit(tT[tT$Gene!="", ]) # remove blanks & NAs
        colnames(tT)[1] <- "Log2FC"
        #write.csv(weights, file="Full_gene_signature.csv", row.names=FALSE)
        filtered.weights.q0.05 <- weights[weights$adj.P.Val<0.05,] # 31 genes
        avg.filtered.weights.q0.05 <- ddply(filtered.weights.q0.05, .(Gene), summarize,
                                            Log2FC = mean(Log2FC, na.rm=TRUE),
                                            sd.Log2FC = sd(Log2FC, na.rm=TRUE)) # prevent duplicate gene names; 31 genes
        
        
        deg.ptm[[types[i]]] <- limma::eBayes()
        deg[[types[i]]] <- limma::eBayes()
      } else {
        deg[[types[i]]] <- limma::eBayes()
      }
    }
    
    #### Step 2. Enrichment analyses ####
    ### ssGSEA ###
    if (!is.null(gmt.genes)) {
      ssGSEA.results <- mGSEA(deg, gmt.genes, FDR = FDR, 
                              min.per.set = min.per.set) 
    }
    
    ### KSEA if relevant ###
    if (any(grepl("phospho", types, ignore.case = TRUE)) & !is.null(gtm.ptm)) {
      KSEA.results <- mGSEA(deg.ptm, gmt.ptm, FDR=FDR, 
                            min.per.set = min.per.set)
    }
    
    ### DMEA ###
    if (!is.null(drug.sensitivity) & !is.null(expression) 
        & !is.null(gmt.drugs)) {
      DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, 
                            expression, deg, FDR = FDR) 
    } 
    
    #### Step 3. Compile results ####
    compiled <- compilEA(deg, ssGSEA.results, DMEA.results, 
                         deg.ptm, KSEA.results)
    
    #### Step 4. Network graphs ####
    networks <- netEA(compiled)
  } else if (length(conditions) == 1) {
    for (j in 2:ncol(data)) {
      #### Step 1. Extract omics data for each sample ####
      temp.data <- list()
      for (i in 1:length(types)) {
        temp.data[[types[i]]] <- data[[types[i]]][ , c(1,j)]
      }
      
      #### Step 2. Enrichment analyses ####
      ### ssGSEA ###
      if (!is.null(gmt.genes)) {
        ssGSEA.results <- mGSEA(temp.data, gmt.genes, FDR = FDR, 
                                min.per.set = min.per.set) 
      }
      
      ### KSEA if relevant ###
      if (any(types %like% "phospho") & !is.null(gmt.ptm)) {
        KSEA.results <- mGSEA(temp.data, gmt.ptm, FDR=FDR, 
                              min.per.set = min.per.set)
      }
      
      #### Step 3. Compile results ####
      compiled <- compilEA(deg, ssGSEA.results, DMEA.results, 
                           deg.ptm, KSEA.results)
      
      #### Step 3. Network graphs ####
      networks <- netEA(compiled)
    }
    ### DMEA ###
    if (!is.null(drug.sensitivity) & !is.null(expression) 
        & !is.null(gmt.drugs)) {
      DMEA.results <- mDMEA(drug.sensitivity, gmt.drugs, expression, 
                            temp.data, FDR = FDR, min.per.set = min.per.set)
    } 
  }
  
  return(list(compiled.results = compiled, network.plots = networks))
}
