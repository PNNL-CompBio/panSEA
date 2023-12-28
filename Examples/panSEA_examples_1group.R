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

##### Step 1: Import CCLE data (PRISM drug AUC, drug moa gmt, RNAseq, proteomics, metabolomics) #####
#### cell line info (CCLE 19Q4)
cell.line.info <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/CCLE_sample_info.csv")
cell.line.info$X <- NULL

#### drug sensitivity (PRISM AUC) for 481 cancer cell lines
AUC.df <- read.csv(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv")
AUC.df$X <- NULL

#### drug sets (PRISM moa)
gmt <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

#### RNAseq (CCLE 19Q4) for 327 adherent cancer cell lines
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
RNA.df <- rbind(RNA.first200, RNA.rest)

#### transcriptomics (IMPROVE)
options(timeout = 360)
download.file("https://figshare.com/ndownloader/files/42131304", "transcriptomics.csv.gz")
transcr.df <- read.csv(gzfile("transcriptomics.csv.gz"),fileEncoding="UTF-16LE")

# merge with genes, samples to stop using improve IDs
transcr.df <- merge(transcr.df, genes)
transcr.df <- merge(transcr.df, samples)

#### proteomics (IMPROVE)
download.file("https://figshare.com/ndownloader/files/41466702", "proteomics.csv.gz")
prot.df <- read.csv(gzfile("proteomics.csv.gz"),fileEncoding="UTF-16LE")

allgenes = readr::read_csv("https://figshare.com/ndownloader/files/40576109")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()

allsamples = readr::read_csv('https://figshare.com/ndownloader/files/40576103')
samples = allsamples|>
  dplyr::select(other_id,improve_sample_id)|>
  unique()

# merge prot.df with genes, samples to stop using improve IDs
prot.df <- merge(prot.df, genes)
prot.df <- merge(prot.df, samples)

#### metabolomics (CCLE)
met.df <- read.csv(paste0("~/OneDrive - PNNL/Documents/GitHub/panSEA/Examples/",
                          "Inputs/CCLE_metabolomics_20190502.csv"))

Sys.setenv("VROOM_CONNECTION_SIZE"=300000) # increase connection buffer size
##### Step 2: Coldren et al (EGFR inhibitor ) #####
### Identify sensitive (SENS) or resistant (RES) cells
### based on Coldren et al, 2006 paper
### not in CCLE RNAseq: CALU3
SENS.cells <- c("NCIH358","NCIH322C","CALU3","NCIH1334","NCIH1648",
                "HCC827","HCC78","NCIH2126","HCC193","HCC95",
                "NCIH3255","HCC4006","NCIH2009","NCIH1650","NCIH820","NCIH1975")
RES.cells <- c("NCIH125","NCIH1703","A549","NCIH157","NCIH460","NCIH520",
               "HCC44","HCC15","NCIH157","NCIH1299","HCC2279")
SENS.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% SENS.cells,] # 5 cell lines: NCIH1650, HCC95, NCIH1975, NCIH1648, NCIH2126
SENS.CCLE_ID <- SENS.cell.info$CCLE_ID
RES.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% RES.cells,] # 7 cell lines: NCIH520, NCIH460, NCIH1299, HCC44, A549, NCIH1703, HCC15
RES.CCLE_ID <- RES.cell.info$CCLE_ID
Coldren.samples <- list(2:6, 7:13)

### Get data for SENS and RES cells
## transcriptomics
SENS.data <- merge(SENS.cell.info, transcr.df, by="CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info, transcr.df, by="CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.full.transcr.data <- cbind(data.SENS,data.RES)

## proteomics
SENS.data <- merge(SENS.cell.info, prot.df, by="CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info, prot.df, by="CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.full.prot.data <- cbind(data.SENS,data.RES)

## metabolomics
SENS.data <- merge(SENS.cell.info, met.df, by="CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info, met.df, by="CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.full.met.data <- cbind(data.SENS,data.RES)

#### Get smaller data sets for comparison to GSE31625
### randomly sample 3 sensitive cell lines & 2 resistant cell lines
sampled.SENS.CCLE_ID <- sample(SENS.CCLE_ID, 3)

sampled.RES.CCLE_ID <- sample(RES.CCLE_ID, 2)

### Get data for SENS and RES cells
## transcriptomics
SENS.data <- merge(SENS.cell.info[SENS.cell.info$CCLE_ID %in% 
                                    sampled.SENS.CCLE_ID, ], transcr.df,
                   by = "CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info[RES.cell.info$CCLE_ID %in% 
                                   sampled.RES.CCLE_ID, ], transcr.df,
                  by = "CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.transcr.data <- cbind(data.SENS,data.RES)

## proteomics
SENS.data <- merge(SENS.cell.info[SENS.cell.info$CCLE_ID %in% 
                                    sampled.SENS.CCLE_ID, ], prot.df,
                   by = "CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info[RES.cell.info$CCLE_ID %in% 
                                  sampled.RES.CCLE_ID, ], prot.df,
                  by = "CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.prot.data <- cbind(data.SENS,data.RES)

## metabolomics
SENS.data <- merge(SENS.cell.info[SENS.cell.info$CCLE_ID %in% 
                                    sampled.SENS.CCLE_ID, ], met.df,
                   by = "CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info[RES.cell.info$CCLE_ID %in% 
                                  sampled.RES.CCLE_ID, ], met.df,
                  by = "CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
Coldren.met.data <- cbind(data.SENS,data.RES)

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

#### Get smaller data sets for comparison to GSE31625
### randomly sample 3 sensitive cell lines & 2 resistant cell lines

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

# get metabolite set info
# based on hmdb v5.0 downloaded 20231220: https://hmdb.ca/downloads 
hmdb <- read.csv(paste0("https://raw.githubusercontent.com/",
                        "BelindaBGarana/panSEA/shiny-app/data/",
                        "ksdb_20231101.csv"))
met.gmt <- DMEA::as_gmt(hmdb, "name", "direct_parent", min.per.set)

#### Step 5a: Fig2: 3 transcriptomic signatures of EGFRi sensitivity ####
types <- c("GSE12790", "GSE31625", "Coldren et al")
data.list <- list(GSE12790, GSE31625, Coldren.transcr)
gmt.features <- list(gene.gmt, gene.gmt, gene.gmt)

Fig2.panSEA <- panSEA::panSEA(data.list, types, group.names = group.names,
                              group.samples = group.samples)

#### Step 5b: Fig3: transcriptomics, proteomics, and metabolomics signatures ####
types <- c("transcriptomics", "proteomics", "metabolomics")
data.list <- list(Coldren.transcr, Coldren.prot, Coldren.met)
feature.names <- c("Gene", "Gene", "Metabolite")
gmt.features <- list(gene.gmt, gene.gmt, met.gmt)

Fig3.panSEA <- panSEA::panSEA(data.list, types, feature.names, 
                              group.names = group.names,
                              group.samples = group.samples,
                              gmt.features = gmt.features)

#### Step 5c: Fig4: Fig3 but only querying non-small cell lung cancer cell line data ####
lung.samples <- cell.line.info[cell.line.info$tissue == "lung", ] # check column, lung filter
lung.RNAseq <- RNAseq[RNAseq$CCLE_ID %in% lung.samples, ]
expression <- list(lung.RNAseq, lung.RNAseq, lung.RNAseq)

Fig3.panSEA <- panSEA::panSEA(data.list, types, feature.names, 
                              group.names = group.names,
                              group.samples = group.samples,
                              gmt.features = gmt.features,
                              expression = expression)
