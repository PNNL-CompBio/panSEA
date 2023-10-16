# Multiomic DMEA
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
# Created: 2023-10-16
# Last modified: 2023-10-16

# input: molecular signature(s) of interest, molecular screening data set(s), drug screen data set, drug annotations
# output: DMEA results for each molecular signature & overall

mDMEA <- function(drug.sensitivity, gmt=NULL, expression, weights, types, value="AUC",
                  sample.names=colnames(expression)[1], gene.names=colnames(weights)[1],
                  weight.values=colnames(weights)[2], rank.metric="Pearson.est", FDR=0.25,
                  num.permutations=1000, stat.type="Weighted", drug.info=NULL, drug="Drug",
                  set.type="moa", min.per.set=6, sep="[|]",
                  exclusions=c("-666", "NA", "na", "NaN", "NULL"), descriptions=NULL,
                  min.per.corr=3, scatter.plots=TRUE, scatter.plot.type="pearson",
                  FDR.scatter.plots=0.05, xlab="Weighted Voting Score", ylab=value,
                  position.x="min", position.y="min", se=TRUE){
  #### Step 1. Check if formats are correct ####
  # check that sample names match between drug.sensitivity and expression data frames
  if (!(sample.names %in% tidyselect::all_of(names(drug.sensitivity),
                                             names(expression)))) {
    stop("sample.names must match across drug.sensitivity and expression data frames")
  }
  
  for (i in 1:length(types)) {
    # check that sample names match 
  }
  
  #### Step 2. Perform DMEA on each omics type ####
  
  #### Step 3. Compile results across omics types ####
}

# filter diffexp for all adjusted p-values < 0.05
contrast.df <- contrast.df[contrast.df$limma_adj < 0.05 & 
                             contrast.df$t_test_adj < 0.05 & 
                             contrast.df$welch_adj < 0.05, ]

## run DMEA wth proteomic BeatAML data
# make sure contrast.df only contains gene symbols in global.df
global.contrast <- contrast.df[contrast.df$feature %in% 
                                 names(global.df)[2:ncol(global.df)], ]

# run DMEA for each contrast
global.contrast.DMEA[[exp.name]] <- DMEA::DMEA(drug.sensitivity = drug.df, 
                                               expression = global.df, weights = global.contrast, 
                                               drug.info = moa.df, sep = ", ")

# merge correlation results with moa annotations for easier query
global.corr <- merge(global.contrast.DMEA[[exp.name]]$corr.result, moa.df)

# save DMEA results locally
setwd(DMEAwd)
dir.create("Proteomic")
setwd("Proteomic")
dir.create(contrasts[i])
setwd(contrasts[i])
dir.create(exp.name)
setwd(exp.name)

# save CSV files
utils::write.csv(global.contrast.DMEA[[exp.name]]$WV.scores, 
                 paste0("DMEA_WV_scores_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(global.contrast.DMEA[[exp.name]]$unused.weights, 
                 paste0("DMEA_unused_weights_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(global.corr, 
                 paste0("DMEA_correlations_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(global.contrast.DMEA[[exp.name]]$result, 
                 paste0("DMEA_results_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(global.contrast.DMEA[[exp.name]]$removed.sets, 
                 paste0("DMEA_removed_sets_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(global.contrast.DMEA[[exp.name]]$unannotated.drugs, 
                 paste0("DMEA_unannotated_drugs_", Sys.Date(), ".csv"), 
                 row.names = FALSE)

# save plots as PDFs
ggplot2::ggsave(paste0("DMEA_volcano_plot_", Sys.Date(), ".pdf"), 
                global.contrast.DMEA[[exp.name]]$volcano.plot, 
                device="pdf")

ggplot2::ggsave(paste0("DMEA_correlation_scatter_plots_", Sys.Date(), ".pdf"), 
                global.contrast.DMEA[[exp.name]]$corr.scatter.plots, device="pdf")

dir.create("mtn_plots")
setwd("mtn_plots")
for(j in 1:length(global.contrast.DMEA[[exp.name]]$mtn.plots)){
  temp.name <- paste0("DMEA_mountain_plot_", 
                      as.filename(names(global.contrast.DMEA[[exp.name]]$mtn.plots)[[j]]), "_", 
                      Sys.Date(), ".pdf")
  ggplot2::ggsave(temp.name, global.contrast.DMEA[[exp.name]]$mtn.plots[[j]], device="pdf")
}

## run DMEA wth proteomic BeatAML data
# make sure contrast.df only contains gene symbols in global.df
RNAseq.contrast <- contrast.df[contrast.df$feature %in% 
                                 names(RNAseq.df)[2:ncol(RNAseq.df)], ]

# run DMEA for each contrast
RNAseq.contrast.DMEA[[exp.name]] <- DMEA::DMEA(drug.sensitivity = drug14.df, 
                                               expression = RNAseq.df, weights = RNAseq.contrast, 
                                               drug.info = moa.df, sep = ", ")

# merge correlation results with moa annotations for easier query
RNAseq.corr <- merge(RNAseq.contrast.DMEA[[exp.name]]$corr.result, moa.df)

# save DMEA results locally
setwd(DMEAwd)
dir.create("Transcriptomic")
setwd("Transcriptomic")
dir.create(contrasts[i])
setwd(contrasts[i])
dir.create(exp.name)
setwd(exp.name)

# save CSV files
utils::write.csv(RNAseq.contrast.DMEA[[exp.name]]$WV.scores, 
                 paste0("DMEA_WV_scores_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(RNAseq.contrast.DMEA[[exp.name]]$unused.weights, 
                 paste0("DMEA_unused_weights_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(RNAseq.corr, 
                 paste0("DMEA_correlations_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(RNAseq.contrast.DMEA[[exp.name]]$result, 
                 paste0("DMEA_results_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(RNAseq.contrast.DMEA[[exp.name]]$removed.sets, 
                 paste0("DMEA_removed_sets_", Sys.Date(), ".csv"), 
                 row.names = FALSE)
utils::write.csv(RNAseq.contrast.DMEA[[exp.name]]$unannotated.drugs, 
                 paste0("DMEA_unannotated_drugs_", Sys.Date(), ".csv"), 
                 row.names = FALSE)

# save plots as PDFs
ggplot2::ggsave(paste0("DMEA_volcano_plot_", Sys.Date(), ".pdf"), 
                RNAseq.contrast.DMEA[[exp.name]]$volcano.plot, 
                device="pdf")

ggplot2::ggsave(paste0("DMEA_correlation_scatter_plots_", Sys.Date(), ".pdf"), 
                RNAseq.contrast.DMEA[[exp.name]]$corr.scatter.plots, device="pdf")

dir.create("mtn_plots")
setwd("mtn_plots")
for(j in 1:length(RNAseq.contrast.DMEA[[exp.name]]$mtn.plots)){
  temp.name <- paste0("DMEA_mountain_plot_", 
                      as.filename(names(RNAseq.contrast.DMEA[[exp.name]]$mtn.plots)[[j]]), "_", 
                      Sys.Date(), ".pdf")
  ggplot2::ggsave(temp.name, RNAseq.contrast.DMEA[[exp.name]]$mtn.plots[[j]], device="pdf")
}  