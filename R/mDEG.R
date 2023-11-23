mDEG <- function(data, types, group.names=c("Diseased", "Healthy"), 
                   group.samples=c(1:(0.5*ncol(data[[1]])), 
                                   (0.5*ncol(data[[1]])+1):ncol(data[[1]]))){
  #### Step 1. Check if formats are correct ####
  # check that there are as many types as data inputs
  if (length(types) != length(data)) {
    stop("Length of types vector must match that of data list")
  }
  
  if (length(conditions) != 2) {
    stop("Only 2 conditions are allowed")
  } else {
    #### Step 2. Differential expression analysis ####
    deg <- list()
    allLogC <- c()
    for (i in 1:length(types)) {
      # make sure there are no duplicated feature names
      if (nrow(data[[i]]) > unique(rownames(data[[i]]))) {
        stop(paste("There must be only one expression value per sample", 
                   "for each feature"))
      } else {
        # separate data based on group (i.e., categorical phenotype)
        data1 <- data[[i]][ , group.samples[1]]
        data2 <- data[[i]][ , group.samples[2]]
        all.data <- cbind(data1, data2)
        
        # create expression sets
        eset <- Biobase::ExpressionSet(all.data)
        
        # log2 transformation if distribution isn't normal
        ex <- exprs(eset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), 
                                  na.rm=TRUE))
        LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(eset) <- log2(ex) }
        allLogC <- c(allLogC, LogC)
        
        # identify sample phenotypes
        gs <- factor(c(rep(0, ncol(data1)), rep(1, ncol(data2))))
        levels(gs) <- make.names(group.names)
        eset$group <- gs
        
        # set up design matrix
        design <- stats::model.matrix(~group + 0, eset)
        colnames(design) <- levels(gs)
        
        # fit linear model
        fit <- limma::lmFit(eset, design)
        
        # set up contrasts of interest and redo fit
        cts <- paste(groups[1], groups[2], sep = "-")
        cont.matrix <- limma::makeContrasts(contrasts=cts, levels=design)
        fit <- limma::contrasts.fit(fit, cont.matrix)
        
        # calculate differential expression and statistics
        fit <- limma::eBayes(fit, 0.01)
        deg[[types[i]]] <- limma::topTable(fit, adjust = "fdr", 
                                           number = nrow(fit))
        deg[[types[i]]]$Gene <- rownames(deg[[types[i]]])
        colnames(deg[[types[i]]])[1] <- "Log2FC"
        deg[[types[i]]] <- 
          deg[[types[i]]][ , c("Gene", 
                               colnames(deg[[types[i]]])
                               [1:(ncol(deg[[types[i]]])-1)])]
      }
    }
  }

  return(list(DEGs = deg,
              Log2Transformed = allLogC))
}
