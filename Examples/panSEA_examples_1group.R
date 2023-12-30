# panSEA examples
# Author: Belinda B. Garana (BG)
# Date created: 2023-12-21
# Last edit: BG 2023-12-21

# To recreate: press run through the whole script
# Overview
# 1 - import CCLE data
# 2 - import Coldren et al. data & run differential expression analysis
# 3 - import GSE31625 data & run differential expression analysis
# 4 - import GSE12790 data & run differential expression analysis
# 5 - run panSEA

rm(list=ls(all=TRUE))
library(plyr);library(dplyr);library(GSA);library(panSEA);library(GEOquery)

# function to save results
savePanSEA <- function(panSEA.results, fig.folder) {
  # set file names
  DMEA.files <- list("DMEA_results.csv" =
                       panSEA.results$mDMEA.results[[1]]$compiled.results$results,
                     "DMEA_mean_results.csv" =
                       panSEA.results$mDMEA.results[[1]]$compiled.results$mean.results,
                     "DMEA_correlation_matrix.pdf" =
                       panSEA.results$mDMEA.results[[1]]$compiled.results$corr.matrix,
                     "DMEA_venn_diagram.pdf" = 
                       panSEA.results$mDMEA.results[[1]]$compiled.results$venn.diagram,
                     "DMEA_dot_plot.pdf" =
                       panSEA.results$mDMEA.results[[1]]$compiled.results$dot.plot,
                     "DMEA_interactive_network.graph.html" =
                       panSEA.results$mDMEA.network[[1]]$interactive)
  GSEA.files <- list("GSEA_results.csv" =
                       panSEA.results$mGSEA.results[[1]]$compiled.results$results,
                     "GSEA_mean_results.csv" =
                       panSEA.results$mGSEA.results[[1]]$compiled.results$mean.results,
                     "GSEA_correlation_matrix.pdf" =
                       panSEA.results$mGSEA.results[[1]]$compiled.results$corr.matrix,
                     "GSEA_venn_diagram.pdf" = 
                       panSEA.results$mGSEA.results[[1]]$compiled.results$venn.diagram,
                     "GSEA_dot_plot.pdf" =
                       panSEA.results$mGSEA.results[[1]]$compiled.results$dot.plot,
                     "GSEA_interactive_network.graph.html" =
                       panSEA.results$mGSEA.network[[1]]$interactive)
  all.files <- list('GSEA' = GSEA.files,
                    'DMEA' = DMEA.files)
  
  base.path <- getwd()
  dir.create(fig.folder)
  setwd(fig.folder)
  for (i in 1:length(all.files)) {
    # create local folder for subset of results
    setwd(paste0(base.path, "/", fig.folder))
    dir.create(names(all.files)[i])
    setwd(names(all.files)[i])
    
    # save results locally
    temp.files <- all.files[[i]]
    
    CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
    for (j in 1:length(CSV.files)) {
      write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
    }
    
    PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
    for (j in 1:length(PDF.files)) {
      ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], device = "pdf")
    }
    
    HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
    if (length(HTML.files) > 0) {
      for (j in 1:length(HTML.files)) {
        visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j])
      }
    }
  }
  setwd(paste0(base.path, "/", fig.folder))
  saveRDS(panSEA.results, paste0("panSEA_", fig.folder, ".rds"))
}

##### Step 1: Import CCLE data (PRISM drug AUC, drug moa gmt, RNAseq, proteomics, metabolomics) #####
#### cell line info (CCLE 19Q4)
cell.line.info <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/CCLE_sample_info.csv")
cell.line.info$X <- NULL

#### drug sensitivity (PRISM AUC) for 481 cancer cell lines
# AUC.df <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv")
# AUC.df$X <- NULL

# #### drug sets (PRISM moa)
# gmt <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

#### RNAseq (CCLE 19Q4) for 327 adherent cancer cell lines
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
RNA.df <- rbind(RNA.first200, RNA.rest)
RNA.first200 <- NULL
RNA.rest <- NULL

#### proteomics (CCLE from IMPROVE)
#download.file("https://figshare.com/ndownloader/files/41466702", "proteomics.csv.gz")
prot.df <- read.csv(gzfile("proteomics.csv.gz"),fileEncoding="UTF-16LE")

allgenes = readr::read_csv("https://figshare.com/ndownloader/files/40576109")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()
genes <- genes[genes$gene_symbol %in% colnames(RNA.df)[2:ncol(RNA.df)], ]

allsamples = readr::read_csv('https://figshare.com/ndownloader/files/40576103')
CCLE.samples <- dplyr::distinct(allsamples[allsamples$id_source == "CCLE",
                                           c("other_id","improve_sample_id")])

# merge prot.df with genes, samples to stop using improve IDs
prot.df <- merge(prot.df, CCLE.samples)
prot.df <- merge(prot.df, genes)
prot.df$entrez_id <- NULL
prot.df <- dplyr::distinct(prot.df)

# convert to wide format for DMEA
prot.df <- reshape2::dcast(prot.df, other_id ~ gene_symbol, mean,
                                value.var = "proteomics")
colnames(prot.df)[1] <- "CCLE_ID"

#### metabolomics (CCLE)
met.df <- read.csv(paste0("~/OneDrive - PNNL/Documents/GitHub/panSEA/Examples/",
                          "Inputs/CCLE_metabolomics_20190502.csv"))
met.df$DepMap_ID <- NULL

Sys.setenv("VROOM_CONNECTION_SIZE"=300000) # increase connection buffer size
##### Step 2: Coldren et al (EGFR inhibitor ) #####
### Identify sensitive (SENS) or resistant (RES) cells
### based on Coldren et al, 2006 paper
SENS.cells <- c("NCIH358","NCIH322C","CALU3","NCIH1334","NCIH1648",
                "HCC827","HCC78","NCIH2126","HCC193","HCC95",
                "NCIH3255","HCC4006","NCIH2009","NCIH1650","NCIH820","NCIH1975")
RES.cells <- c("NCIH125","NCIH1703","A549","NCIH157","NCIH460","NCIH520",
               "HCC44","HCC15","NCIH157","NCIH1299","HCC2279")
SENS.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% SENS.cells,] 
SENS.stripped_cell_line_name <- SENS.cell.info$stripped_cell_line_name
# 12 SENS cell lines: "HCC827", "NCIH1650", "HCC4006", "NCIH3255", "CALU3",
# "HCC95", "HCC78", "NCIH1975", "NCIH1648", "NCIH2126", "NCIH358", "NCIH2009"

RES.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% RES.cells,]
RES.stripped_cell_line_name <- RES.cell.info$stripped_cell_line_name
# 8 RES cell lines: "NCIH520", "NCIH460", "NCIH1299", "HCC44", "A549", 
# "HCC2279", "NCIH1703", "HCC15"

# Coldren.samples <- list(seq(2, length.out = length(SENS.stripped_cell_line_name)), 
#                         seq((length(SENS.stripped_cell_line_name) + 2), 
#                             length.out = length(RES.stripped_cell_line_name)))

Coldren.sample.IDs <- CCLE.samples[CCLE.samples$other_id %in% 
                                     c(SENS.cell.info$CCLE_Name,
                                       RES.cell.info$CCLE_Name),]
colnames(Coldren.sample.IDs)[1] <- "CCLE_Name" # change other_id to CCLE_Name
Coldren.sample.IDs <- merge(Coldren.sample.IDs, cell.line.info)
# Coldren.SENS.IDs <- Coldren.sample.IDs[Coldren.sample.IDs$CCLE_Name %in% 
#                                          SENS.cell.info$CCLE_Name, 
#                                        c("improve_sample_id", "CCLE_ID")]
# Coldren.RES.IDs <- Coldren.sample.IDs[Coldren.sample.IDs$CCLE_Name %in% 
#                                         RES.cell.info$CCLE_Name, 
#                                       c("improve_sample_id", "CCLE_ID")]

### Get data for SENS and RES cells
## transcriptomics
SENS.data <- RNA.df[RNA.df$CCLE_ID %in% SENS.cell.info$CCLE_Name, ] 
# 5 SENS cell lines: "NCIH1650_LUNG", "HCC95_LUNG", "NCIH1975_LUNG",
# "NCIH1648_LUNG", "NCIH2126_LUNG"
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL

RES.data <- RNA.df[RNA.df$CCLE_ID %in% RES.cell.info$CCLE_Name, ] 
# 7 cell lines: "NCIH520_LUNG", "NCIH460_LUNG", "NCIH1299_LUNG", "HCC44_LUNG",
# "A549_LUNG", "NCIH1703_LUNG", "HCC15_LUNG" 
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL

data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.transcr <- as.data.frame(cbind(data.SENS,data.RES))
Coldren.transcr$Gene <- rownames(Coldren.transcr)
Coldren.transcr <- Coldren.transcr[ , c("Gene", 
                                        colnames(Coldren.transcr)[
                                          1:(ncol(Coldren.transcr)-1)])]

## proteomics
SENS.data <- prot.df[prot.df$CCLE_ID %in% SENS.cell.info$CCLE_Name, ] 
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
# 8 cell lines: "HCC827_LUNG", "HCC95_LUNG", "NCIH1650_LUNG", "NCIH1975_LUNG",
# "NCIH2009_LUNG", "NCIH2126_LUNG", "NCIH3255_LUNG", "NCIH358_LUNG" 

RES.data <- prot.df[prot.df$CCLE_ID %in% RES.cell.info$CCLE_Name, ] 
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
# 7 cell lines: "A549_LUNG", "HCC15_LUNG", "HCC44_LUNG", "NCIH1299_LUNG",
# "NCIH1703_LUNG", "NCIH460_LUNG", "NCIH520_LUNG" 

data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.prot <- as.data.frame(cbind(data.SENS,data.RES))
Coldren.prot$Gene <- rownames(Coldren.prot)
Coldren.prot <- Coldren.prot[ , c("Gene", 
                                  colnames(Coldren.prot)[
                                    1:(ncol(Coldren.prot)-1)])]

## metabolomics
SENS.data <- met.df[met.df$CCLE_ID %in% SENS.cell.info$CCLE_Name, ] 
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
# 12 cell lines: "NCIH358_LUNG", "NCIH1650_LUNG", "HCC95_LUNG",
# "NCIH3255_LUNG", "NCIH1648_LUNG", "NCIH1975_LUNG", "HCC78_LUNG",
# "NCIH2126_LUNG", "HCC827_LUNG", "CALU3_LUNG", "NCIH2009_LUNG", "HCC4006_LUNG" 

RES.data <- met.df[met.df$CCLE_ID %in% RES.cell.info$CCLE_Name, ] 
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
# 8 cell lines: "HCC2279_LUNG", "NCIH520_LUNG", "HCC15_LUNG", "NCIH460_LUNG",
# "A549_LUNG", "HCC44_LUNG", "NCIH1703_LUNG", "NCIH1299_LUNG"

data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.met <- as.data.frame(cbind(data.SENS,data.RES))
Coldren.met$Metabolite <- rownames(Coldren.met)
Coldren.met <- Coldren.met[ , c("Metabolite", 
                                colnames(Coldren.met)[
                                  1:(ncol(Coldren.met)-1)])]

### run differential expression analysis
### have to do this separately for each omics because they have different
### numbers of SENS and RES samples
## transcriptomics
# group.samples <- list(seq(2, length.out = nrow(RNA.df[
#   RNA.df$CCLE_ID %in% SENS.cell.info$CCLE_Name, ])), 
#   seq((nrow(RNA.df[RNA.df$CCLE_ID %in% 
#                          SENS.cell.info$CCLE_Name, ]) + 2), 
#       length.out = nrow(met.df[met.df$CCLE_ID %in% RES.cell.info$CCLE_Name, ])))
group.samples <- list(RNA.df$CCLE_ID[RNA.df$CCLE_ID %in% 
                                       SENS.cell.info$CCLE_Name],
                      RNA.df$CCLE_ID[RNA.df$CCLE_ID %in% 
                                       RES.cell.info$CCLE_Name])
Coldren.transcr.DEGs <- panSEA::mDEG(list(Coldren.transcr), types = "Transcriptomics",
                                     group.names = c("Sensitive", "Resistant"),
                                     group.samples, "Gene")$all.results[[1]]
write.csv(Coldren.transcr.DEGs, "Coldren_transcriptomics_DEGs_before_avg.csv", 
          row.names = FALSE)
Coldren.transcr.DEGs <- read.csv("Coldren_transcriptomics_DEGs_before_avg.csv")
library(plyr);library(dplyr)
Coldren.transcr <- ddply(Coldren.transcr.DEGs, .(Gene), summarise,
                                    Log2FC = mean(Log2FC, na.rm=TRUE)) #prevent duplicate gene names
write.csv(Coldren.transcr, "Coldren_transcriptomics_DEGs.csv", 
          row.names = FALSE)
Coldren.transcr <- read.csv("Coldren_transcriptomics_DEGs.csv")

## proteomics
group.samples <- list(prot.df$CCLE_ID[prot.df$CCLE_ID %in% 
                                       SENS.cell.info$CCLE_Name],
                      prot.df$CCLE_ID[prot.df$CCLE_ID %in% 
                                       RES.cell.info$CCLE_Name])
Coldren.prot.DEGs <- panSEA::mDEG(list(Coldren.prot), types = "Proteomics",
                                     group.names = c("Sensitive", "Resistant"),
                                     group.samples, "Gene")$all.results[[1]]
Coldren.prot.DEGs <- plyr::ddply(Coldren.prot.DEGs, .(Gene), summarise,
                                     Log2FC = mean(Log2FC, na.rm=TRUE)) #prevent duplicate gene names
write.csv(Coldren.prot.DEGs, "Coldren_proteomics_DEGs.csv", 
          row.names = FALSE)
Coldren.prot.DEGs <- read.csv("Coldren_proteomics_DEGs.csv")

## metabolomics
group.samples <- list(met.df$CCLE_ID[met.df$CCLE_ID %in% 
                                       SENS.cell.info$CCLE_Name],
                      met.df$CCLE_ID[met.df$CCLE_ID %in% 
                                       RES.cell.info$CCLE_Name])
Coldren.met.DEGs <- panSEA::mDEG(list(Coldren.met), types = "Metabolomics",
                                     group.names = c("Sensitive", "Resistant"),
                                     group.samples, "Metabolite")$all.results[[1]]
Coldren.met.DEGs <- plyr::ddply(Coldren.met.DEGs, .(Metabolite), summarise,
                                 Log2FC = mean(Log2FC, na.rm=TRUE)) #prevent duplicate gene names
write.csv(Coldren.met.DEGs, "Coldren_metabolomics_DEGs.csv", 
          row.names = FALSE)
Coldren.met.DEGs <- read.csv("Coldren_metabolomics_DEGs.csv")

##### Step 3: GSE31625 (EGFR inhibitor) #####
# load series and platform data from GEO
gset <- getGEO("GSE31625", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("X11111111XXXXXXXXX0000000XXXXXXXX0000X0000001111") # 0 is phenotype of numerator in FC; 1 of denominator; X: excluded
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

### run differential expression analysis
# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
SENS.expr <- as.data.frame(exprs(gset[,which(sml==0)]))
RES.expr <- as.data.frame(exprs(gset[,which(sml==1)]))

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("test","test2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- limma::lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- limma::makeContrasts(contrasts=cts, levels=design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- limma::eBayes(fit2, 0.01)
tT <- limma::topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")) # 22,215 genes
tT.noblanks <- tT[tT$Gene.symbol!="",] # 21,146 genes
colnames(tT.noblanks)[6] <- "Log2FC"
colnames(tT.noblanks)[7] <- "Gene"
GSE31625.DEGs <- tT.noblanks
GSE31625.DEGs$Log2Transformed <- LogC
write.csv(GSE31625.DEGs, "GSE31625_transcriptomics_DEGs_before_avg.csv", 
          row.names = FALSE)
GSE31625.DEGs <- read.csv("GSE31625_transcriptomics_DEGs_before_avg.csv")
GSE31625 <- plyr::ddply(GSE31625.DEGs, .(Gene), summarise,
                                 Log2FC = mean(Log2FC, na.rm=TRUE)) #prevent duplicate gene names
write.csv(GSE31625, "GSE31625_transcriptomics_DEGs.csv", 
          row.names = FALSE)
GSE31625 <- read.csv("GSE31625_transcriptomics_DEGs.csv")

# GSE31625.DEGs <- panSEA::mDEG(GSE31625, types = "Transcriptomics", 
#                              group.names = c("Sensitive", "Resistant"),
#                              list(sml[which(sml == 0)], sml[which(sml == 1)])
#                              )$all.results

##### Step 4: GSE12790 (EGFR inhibitor) #####
# load series and platform data from GEO
gset <- getGEO("GSE12790", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX1XX",
               "X0X11X1XX1X1X1XX1XX1X1110X1XX1XXX11X1XX0XX1XXXXX") # 0 is phenotype of numerator in FC; 1 of denominator; X: excluded
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

### run differential expression analysis
# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("test","test2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- limma::lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- limma::makeContrasts(contrasts=cts, levels=design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- limma::eBayes(fit2, 0.01)
tT <- limma::topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title")) # 54,675 genes
tT.noblanks <- tT[tT$Gene.symbol!="",] #45,118 genes
colnames(tT.noblanks)[6] <- "Log2FC"
colnames(tT.noblanks)[7] <- "Gene"
GSE12790.DEGs <- tT.noblanks
GSE12790.DEGs$Log2Transformed <- LogC
write.csv(GSE12790.DEGs, "GSE12790_transcriptomics_DEGs_before_avg.csv", 
          row.names = FALSE)
GSE12790.DEGs <- read.csv("GSE12790_transcriptomics_DEGs_before_avg.csv")
GSE12790 <- dplyr::distinct(plyr::ddply(GSE12790.DEGs, .(Gene), summarise,
                             Log2FC = mean(Log2FC, na.rm=TRUE))) #prevent duplicate gene names
write.csv(GSE12790, "GSE12790_transcriptomics_DEGs.csv", 
          row.names = FALSE)
GSE12790 <- read.csv("GSE12790_transcriptomics_DEGs.csv")

# GSE12790.DEGs <- panSEA::mDEG(GSE12790, types = "Transcriptomics", 
#                               group.names = c("Sensitive", "Resistant"),
#                               list(sml[which(sml == 0)], sml[which(sml == 1)])
#                               )$all.results

##### Step 5: run panSEA #####
# IMPORTANT: each group must have the same # & position of samples across data sets
# e.g. "Sensitive" samples must be in the same columns across input data.list
# but "Resistant" could have different # as "Sensitive" hypothetically
group.names <- "Sensitive vs. Resistant"
group.samples <- 2

# get gene set info
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C2", "CP:KEGG")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gene.gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description")
saveRDS(gene.gmt, "gmt_MSigDB_Homo-sapiens_C2_CP:KEGG.rds")

# get metabolite set info
# based on hmdb v5.0 downloaded 20231220: https://hmdb.ca/downloads 
hmdb <- XML::xmlToDataFrame("hmdb_metabolites.xml")
saveRDS(hmdb, "hmdb_metabolites_dataframe.rds")
hmdb <- readRDS("hmdb_metabolites_dataframe.rds")
#hmdb$class <- sub(".*belongs to the class of organic compounds known as (.+)\\..*", "\\1", hmdb$taxonomy)
hmdb$class <- NA
hmdb$class <- stringr::str_match(hmdb$taxonomy, "belongs to the class of organic compounds known as\\s*(.*?)\\s*\\. ")[ , 2]
#inorganic.hmdb <- hmdb[is.na(hmdb$class), ] # 72291 / 217920
hmdb[is.na(hmdb$class), ]$class <- 
  stringr::str_match(
    hmdb[is.na(hmdb$class), ]$taxonomy, 
    "belongs to the class of inorganic compounds known as\\s*(.*?)\\s*\\. ")[ , 2]
#unclassified.hmdb <- hmdb[is.na(hmdb$class), ] # 71984 / 217920
#nonblank.unclassified.hmdb <- unclassified.hmdb[unclassified.hmdb$taxonomy != "\n  ", ] # 1 / 71984
hmdb[is.na(hmdb$class), ]$class <- 
  stringr::str_match(
    hmdb[is.na(hmdb$class), ]$taxonomy, 
    "belongs to the class of chemical entities known as\\s*(.*?)\\s*\\. ")[ , 2]
#unclassified.hmdb <- hmdb[is.na(hmdb$class), ] # 71983 / 217920
saveRDS(hmdb, "hmdb_metabolites_dataframe_with_class.rds")
hmdb <- dplyr::distinct(na.omit(hmdb[ , c("name", "class")])) # 145937 metabolites

met.gmt <- DMEA::as_gmt(hmdb, "name", "class", min.per.set = 6)
saveRDS(met.gmt, "gmt_HMDB.rds")

hmdb.CCLE <- hmdb[tolower(hmdb$name) %in% colnames(met.df)[2:ncol(met.df)], ] # 46 / 225
hmdb.CCLE$name <- lapply(hmdb.CCLE$name, tolower)
met.gmt <- DMEA::as_gmt(hmdb.CCLE, "name", "class", min.per.set = 6)
saveRDS(met.gmt, "gmt_HMDB_in_CCLE.rds")

#### Step 5a: Fig2: 3 transcriptomic signatures of EGFRi sensitivity ####
types <- c("GSE12790", "GSE31625", "Coldren et al")
data.list <- list("GSE12790" = GSE12790, 
                  "GSE31625" = GSE31625, 
                  "Coldren et al" = Coldren.transcr)
gmt.features <- list(gene.gmt, gene.gmt, gene.gmt)
expression <- list(RNA.df, RNA.df, RNA.df)

Fig2.panSEA <- panSEA::panSEA(data.list, types, group.names = group.names,
                              group.samples = group.samples, gmt.features = gmt.features,
                              expression = expression)
# save results
fig.folder <- "Fig2"
savePanSEA(Fig2.panSEA, fig.folder)
Fig2.panSEA <- readRDS(paste0("panSEA_", fig.folder, ".rds"))
Fig2.panSEA <- NULL # make space to process next analysis

#### Step 5b: Fig3: transcriptomics, proteomics, and metabolomics signatures ####
types <- c("transcriptomics", "proteomics", "metabolomics")
data.list <- list(Coldren.transcr.DEGs, Coldren.prot.DEGs, Coldren.met.DEGs)
feature.names <- c("Gene", "Gene", "Metabolite")
gmt.features <- list(gene.gmt, gene.gmt, met.gmt)
expression <- list(RNA.df, prot.df, met.df)

Fig3.panSEA <- panSEA::panSEA(data.list, types, feature.names, 
                              group.names = group.names,
                              group.samples = group.samples,
                              gmt.features = gmt.features)

# save results
savePanSEA(Fig3.panSEA, "Fig3")
Fig3.panSEA <- NULL # make space to process next analysis

#### Step 5c: Fig4: Fig3 but only querying non-small cell lung cancer cell line data ####
NSCLC.samples <- cell.line.info[cell.line.info$lineage_subtype == "NSCLC", ] # 158 CCLE
NSCLC.RNA.df <- RNA.df[RNA.df$CCLE_ID %in% NSCLC.samples$CCLE_Name, ] # 51 CCLE
NSCLC.prot.df <- prot.df[prot.df$CCLE_ID %in% NSCLC.samples$CCLE_Name, ] # 51 CCLE
NSCLC.met.df <- met.df[met.df$CCLE_ID %in% NSCLC.samples$CCLE_Name, ] # 51 CCLE
expression <- list(NSCLC.RNA.df, NSCLC.prot.df, NSCLC.met.df)

Fig4.panSEA <- panSEA::panSEA(data.list, types, feature.names, 
                              group.names = group.names,
                              group.samples = group.samples,
                              gmt.features = gmt.features,
                              expression = expression)

# save results
savePanSEA(Fig4.panSEA, "Fig4")
Fig4.panSEA <- NULL # make space to process next analysis
