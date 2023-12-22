# panSEA examples
# Author: Belinda B. Garana (BG)
# Date created: 2023-12-21
# Last edit: BG 2023-12-21

# To recreate: press run through the whole script
# Overview
# 1 - import CCLE data
# 2 - import Coldren et al. data
# 3 - import GSE31625 data
# 4 - import GSE12790 data
# 5 - run panSEA

rm(list=ls(all=TRUE))
Sys.setenv("VROOM_CONNECTION_SIZE"=300000) # increase connection buffer size
library(plyr);library(dplyr);library(GSA);library(panSEA);library(GEOquery)

##### Step 1: Import CCLE data (PRISM drug AUC, drug moa gmt, RNAseq, proteomics, metabolomics) #####
#### cell line info (CCLE 19Q4)
cell.line.info <- read.csv(file="https://raw.github.com/BelindaBGarana/panSEA/shiny-app/Inputs/CCLE_sample_info.csv")
cell.line.info$X <- NULL

#### drug sensitivity (PRISM AUC) for 481 cancer cell lines
AUC.df <- read.csv(file="https://raw.github.com/BelindaBGarana/panSEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv")
AUC.df$X <- NULL

#### drug sets (PRISM moa)
gmt <- GSA.read.gmt(file="https://raw.github.com/BelindaBGarana/panSEA/shiny-app/Inputs/MOA_gmt_file_n6_no_special_chars.gmt")

#### RNAseq (CCLE 19Q4) for 327 adherent cancer cell lines
download.file("https://raw.github.com/BelindaBGarana/panSEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file("https://raw.github.com/BelindaBGarana/panSEA/shiny-app/Inputs/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin", 
              destfile = paste0(getwd(),"/Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"))
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
RNA.df <- rbind(RNA.first200, RNA.rest)

#### proteomics (CCLE)

#### metabolomics (CCLE)

##### Step 2: Coldren et al (EGFR inhibitor ) #####
### Identify sensitive (SENS) or resistant (RES) cells
### based on Coldren et al, 2006 paper
### not in CCLE RNAseq: CAL
SENS.cells <- c("NCIH358","NCIH322C","CALU3","NCIH1334","NCIH1648",
                "HCC827","HCC78","NCIH2126","HCC193","HCC95",
                "NCIH3255","HCC4006","NCIH2009","NCIH1650","NCIH820","NCIH1975") # not in CCLE RNAseq: CALU3
RES.cells <- c("NCIH125","NCIH1703","A549","NCIH157","NCIH460","NCIH520",
               "HCC44","HCC15","NCIH157","NCIH1299","HCC2279")
SENS.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% SENS.cells,] # 5 cell lines: NCIH1650, HCC95, NCIH1975, NCIH1648, NCIH2126
SENS.CCLE_ID <- SENS.cell.info$CCLE_ID
RES.cell.info <- cell.line.info[cell.line.info$stripped_cell_line_name %in% RES.cells,] # 7 cell lines: NCIH520, NCIH460, NCIH1299, HCC44, A549, NCIH1703, HCC15
RES.CCLE_ID <- RES.cell.info$CCLE_ID

#### Get data for SENS and RES cells
### transcriptomics
SENS.data <- merge(SENS.cell.info, RNA.df, by="CCLE_ID")
SENS.data$stripped_cell_line_name <- NULL
rownames(SENS.data) <- SENS.data$CCLE_ID
SENS.data$CCLE_ID <- NULL
RES.data <- merge(RES.cell.info, RNA.df, by="CCLE_ID")
RES.data$stripped_cell_line_name <- NULL
rownames(RES.data) <- RES.data$CCLE_ID
RES.data$CCLE_ID <- NULL
data.SENS <- t(SENS.data)
data.RES <- t(RES.data)
SENS.RES.data <- cbind(data.SENS,data.RES)
sml <- c(0,0,0,0,0,1,1,1,1,1,1,1)

### proteomics

### metabolomics

#### Get smaller data sets for comparison to GSE31625
### randomly sample 3 sensitive cell lines & 2 resistant cell lines


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
group.names = c("Sensitive", "Resistant")
#### Step 5a: Fig2: 3 transcriptomic signatures of EGFRi sensitivity ####
types <- c("GSE12790", "GSE31625", "Coldren et al")
data.list <- list()
group.samples <- list()

Fig2.panSEA <- panSEA::panSEA(data.list, types, group.names = group.names,
                              group.samples = group.samples)

#### Step 5b: Fig3: transcriptomics, proteomics, and metabolomics signatures ####
types <- c("transcriptomics", "proteomics", "metabolomics")
data.list <- list()
feature.names <- c("Gene", "Gene", "Metabolite")
group.samples <- list()
gmt.features <- list("msigdb_Homo sapiens_C2_CP:KEGG",
                     "msigdb_Homo sapiens_C2_CP:KEGG",
                     "hmdb")

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
