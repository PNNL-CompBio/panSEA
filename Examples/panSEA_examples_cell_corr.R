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
base.path <- "~/OneDrive - PNNL/Documents/GitHub/panSEA/Examples/"
setwd(paste0(base.path, "Inputs"))

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
download.file("https://figshare.com/ndownloader/files/41466702", "proteomics.csv.gz")
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
setwd(base.path)
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

# compile data into list for panSEA input
all.Coldren.CCLE_IDs <- colnames(Coldren.transcr)[2:ncol(Coldren.transcr)]
all.Coldren.names <- c("H1650", "HCC95", "H1975", "H1648", "H2126", "H520", 
                       "H460", "H1299", "HCC44", "A549", "H1703", "HCC15")
Coldren <- list()
for (i in 2:ncol(Coldren.transcr)) {
  temp.df <- Coldren.transcr[ , c(1, i)]
  colnames(temp.df)[2] <- "Input"
  Coldren[[all.Coldren.names[i-1]]] <- temp.df
}

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

all.Coldren.prot.CCLE_IDs <- colnames(Coldren.prot)[2:ncol(Coldren.prot)]
all.Coldren.prot.names <- c("HCC827", "HCC95", "H1650", "H1975", "H2009",
                            "H2126", "H3255", "H358", "A549", "HCC15", "HCC44",
                            "H1299", "H1703", "H460", "H520")
# compile data into list for panSEA input
Coldren.prot.list <- list()
for (i in 2:ncol(Coldren.prot)) {
  temp.df <- Coldren.prot[ , c(1, i)]
  colnames(temp.df)[2] <- "Input"
  Coldren.prot.list[[all.Coldren.prot.names[i-1]]] <- temp.df
}

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

all.Coldren.met.CCLE_IDs <- colnames(Coldren.met)[2:ncol(Coldren.met)]
all.Coldren.met.names <- c("H358", "H1650", "HCC95", "H3255", "H1648",
                           "H1975", "HCC78", "H2126", "HCC827", "CALU3", 
                           "H2009", "HCC4006", "HCC2279", "H520", "HCC15", 
                           "H460", "A549", "HCC44", "H1703", "H1299")
# compile data into list for panSEA input
Coldren.met.list <- list()
for (i in 2:ncol(Coldren.met)) {
  temp.df <- Coldren.met[ , c(1, i)]
  colnames(temp.df)[2] <- "Input"
  Coldren.met.list[[all.Coldren.met.names[i-1]]] <- temp.df
}

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

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

SENS.expr <- as.data.frame(exprs(gset[,which(sml==0)]))
RES.expr <- as.data.frame(exprs(gset[,which(sml==1)]))
SENS.GSE31625 <- colnames(SENS.expr)
# matches with "H1650-8", "H1650-9", "H1650-H", "H1650A", "H1650B", "H1650C",
# "H1650p8", "H3255_1", "H3255_2", "H3255_3", "H32551-Dec", "PC-9_7", "PC-9_8",
# "PC-9_9", "PC-9_10", "PC9", "PC92-4val" on GEO
SENS.GSE31625.names <- c("H1650_1", "H1650_2", "H1650_3", "H1650_4", "H1650_5", 
                         "H1650_6", "H1650_7", "H3255_1", "H3255_2", "H3255_3", 
                         "H3255_4", "PC9_1", "PC9_2", "PC9_3", "PC9_4", "PC9_5",
                         "PC9_6")
colnames(SENS.expr) <- SENS.GSE31625.names

RES.GSE31625 <- colnames(RES.expr)
# matches with "A549B", "A549C", "A549-3.05", "A549p2", "A549p3", "A549p4", 
# "A549p5", "A549p6", "UKY-29_4", "UKY-29_5", "UKY-29_6", "UKY29-3" on GEO
RES.GSE31625.names <- c("A549_1", "A549_2", "A549-3", "A549_4", "A549_5", 
                        "A549_6", "A549_7", "A549_8", "UKY29_1", "UKY29_2", 
                        "UKY29_3", "UKY29_4")
colnames(RES.expr) <- RES.GSE31625.names

GSE31625.df <- cbind(SENS.expr, RES.expr)
GSE31625.df$Gene <- rownames(GSE31625.df)
GSE31625.df <- GSE31625.df[ , c("Gene", 
                                colnames(GSE31625.df)[
                                  1:(ncol(GSE31625.df)-1)])]

all.GSE31625.names <- c(SENS.GSE31625.names, RES.GSE31625.names)
GSE31625.names <- c("H1650", "H3255", "PC9", "A549", "UKY29")
GSE31625.CCLE_IDs <- c("NCIH1650_LUNG", "NCIH3255_LUNG", "PC9_LUNG", 
                       "A549_LUNG")
all.GSE31625.CCLE_IDs <- c()
for (i in 1:length(all.GSE31625.names)) {
  if (grepl(GSE31625.names[1], all.GSE31625.names[i])) {
    all.GSE31625.CCLE_IDs <- c(all.GSE31625.CCLE_IDs, GSE31625.CCLE_IDs[1])
  } else if (grepl(GSE31625.names[2], all.GSE31625.names[i])) {
    all.GSE31625.CCLE_IDs <- c(all.GSE31625.CCLE_IDs, GSE31625.CCLE_IDs[2])
  } else if (grepl(GSE31625.names[3], all.GSE31625.names[i])) {
    all.GSE31625.CCLE_IDs <- c(all.GSE31625.CCLE_IDs, GSE31625.CCLE_IDs[3])
  } else if (grepl(GSE31625.names[4], all.GSE31625.names[i])) {
    all.GSE31625.CCLE_IDs <- c(all.GSE31625.CCLE_IDs, GSE31625.CCLE_IDs[4])
  } else {
    all.GSE31625.CCLE_IDs <- c(all.GSE31625.CCLE_IDs, NA)
  }
}

# compile data into list for panSEA input
GSE31625 <- list()
for (i in 2:ncol(GSE31625.df)) {
  temp.df <- GSE31625.df[ , c(1, i)]
  colnames(temp.df)[2] <- "Input"
  GSE31625[[colnames(GSE31625.df)[i]]] <- temp.df
}
temp.df <- NULL

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

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
SENS.expr <- as.data.frame(exprs(gset[,which(sml==0)]))
RES.expr <- as.data.frame(exprs(gset[,which(sml==1)]))
SENS.GSE12790 <- colnames(SENS.expr)
# matches with "breast cancer cell line_HCC1806", 
# "breast cancer cell line_CAL85-1", "breast cancer cell line_HDQ-P1" on GEO
SENS.GSE12790.names <- c("HCC1806", "CAL851", "HDQ-P1")
colnames(SENS.expr) <- SENS.GSE12790.names

RES.GSE12790 <- colnames(RES.expr)
# matches with "breast cancer cell line_BT20", "breast cancer cell line_BT474",
# "breast cancer cell line_BT483", "breast cancer cell line_CAMA-1",
# "breast cancer cell line_MCF7", "breast cancer cell line_MDA-MB-231",
# "breast cancer cell line_MDA-MB-231", "breast cancer cell line_MDA-MB-453",
# "breast cancer cell line_T47D", "breast cancer cell line_ZR75-1",
# "breast cancer cell line_CAL-120", "breast cancer cell line_CAL-51", 
# "breast cancer cell line_EFM-192A", "breast cancer cell line_KPL1", 
# "breast cancer cell line_HCC1428", "breast cancer cell line_HCC1954", 
# "breast cancer cell line_HCC1500", "breast cancer cell line_HCC1569" on GEO
RES.GSE12790.names <- c("BT20", "BT474", "BT483", "CAMA-1", "MCF7", 
                        "MDAMB231_1", "MDAMB231_2", "MDAMB453", "T47D", 
                        "ZR751", "CAL120", "CAL51", "EFM192A", "KPL1", 
                        "HCC1428", "HCC1954", "HCC1500", "HCC1569")
colnames(RES.expr) <- RES.GSE12790.names

GSE12790.df <- cbind(SENS.expr, RES.expr)
GSE12790.df$Gene <- rownames(GSE12790.df)
GSE12790.df <- GSE12790.df[ , c("Gene", 
                                colnames(GSE12790.df)[
                                  1:(ncol(GSE12790.df)-1)])]

all.GSE12790.names <- c(SENS.GSE12790.names, RES.GSE12790.names)
GSE12790.names <- c("HCC1806", "CAL851", "HDQ-P1", "BT20", "BT474", "BT483",
                    "CAMA1", "MCF7", "MDAMB231", "MDAMB453", "T47D", 
                    "ZR751", "CAL120", "CAL51", "EFM192A", "KPL1",
                    "HCC1428", "HCC1954", "HCC1500", "HCC1569")
GSE12790.CCLE_IDs <- c("HCC1806_BREAST", "CAL851_BREAST", "HDQP1_BREAST", 
                       "BT20_BREAST", "BT474_BREAST", "BT483_BREAST", 
                       "CAMA1_BREAST", "MCF7_BREAST", "MDAMB231_BREAST",
                       "MDAMB453_BREAST", "T47D_BREAST", "ZR751_BREAST",
                       "CAL120_BREAST", "CAL51_BREAST", "EFM192A_BREAST",
                       "KPL1_BREAST", "HCC1428_BREAST", "HCC1954_BREAST",
                       "HCC1500_BREAST", "HCC1569_BREAST")
all.GSE12790.CCLE_IDs <- c()
for (i in 1:length(all.GSE12790.names)) {
  match <- FALSE
  for (j in 1:length(GSE12790.names)) {
    if (grepl(GSE12790.names[j], all.GSE12790.names[i])) {
      all.GSE12790.CCLE_IDs <- c(all.GSE12790.CCLE_IDs, GSE12790.CCLE_IDs[j])
      match <- TRUE
      break
    }
  }
  if (!match) {
    all.GSE12790.CCLE_IDs <- c(all.GSE12790.CCLE_IDs, NA)
  }
}

# compile data into list for panSEA input
GSE12790 <- list()
for (i in 2:ncol(GSE12790.df)) {
  temp.df <- GSE12790.df[ , c(1, i)]
  colnames(temp.df)[2] <- "Input"
  GSE12790[[colnames(GSE12790.df)[i]]] <- temp.df
}
temp.df <- NULL

##### Step 5: run panSEA #####
dir.create("cell_corr")
setwd("cell_corr")
# IMPORTANT: each group must have the same # & position of samples across data sets
# e.g. "Sensitive" samples must be in the same columns across input data.list
# but "Resistant" could have different # as "Sensitive" hypothetically
group.names <- "Sensitive or Resistant"
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
gene.gmt <- readRDS("Inputs/gmt_MSigDB_Homo-sapiens_C2_CP:KEGG.rds")

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
hmdb <- readRDS("hmdb_metabolites_dataframe_with_class.rds")
hmdb <- dplyr::distinct(na.omit(hmdb[ , c("name", "class")])) # 145937 metabolites

met.gmt <- DMEA::as_gmt(hmdb, "name", "class", min.per.set = 6)
saveRDS(met.gmt, "gmt_HMDB.rds")

hmdb <- readRDS("hmdb_metabolites_dataframe_with_class.rds")
hmdb <- dplyr::distinct(na.omit(hmdb[ , c("name", "class")])) # 145937 metabolites
# hmdb.wide <- reshape2::dcast(hmdb, class ~ name)
# hmdb.match <- colnames(hmdb.wide)[tolower(colnames(hmdb.wide)) %in% colnames(met.df)[2:ncol(met.df)]] # 46 / 225
hmdb.CCLE <- hmdb[tolower(hmdb$name) %in% colnames(met.df)[2:ncol(met.df)], ] # 46 / 225
hmdb.CCLE$name <- lapply(hmdb.CCLE$name, tolower)
met.gmt <- DMEA::as_gmt(hmdb.CCLE, "name", "class", min.per.set = 6)
saveRDS(met.gmt, "gmt_HMDB_in_CCLE.rds")
hmdb <- NULL
met.gmt <- readRDS("gmt_HMDB_in_CCLE.rds")

# increase # of groups with 6+ metabolites
hmdb.CCLE[grep("alpha", hmdb.CCLE$class), ]$class <- "alpha amino acids"
met.gmt <- DMEA::as_gmt(hmdb.CCLE, "name", "class", min.per.set = 6)
saveRDS(met.gmt, "gmt_HMDB_in_CCLE_2sets.rds")
met.gmt <- readRDS("Inputs/gmt_HMDB_in_CCLE_2sets.rds")

# met.gmt <- DMEA::as_gmt(hmdb.CCLE, "name", "class", min.per.set = 4)
# saveRDS(met.gmt, "gmt_HMDB_in_CCLE_4minperset.rds")

#### Step 5a: GSE12790 ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.GSE12790.names)) {
  gmt.features[[all.GSE12790.names[i]]] <- gene.gmt
  expression[[all.GSE12790.names[i]]] <- RNA.df[
    RNA.df$CCLE_ID != all.GSE12790.CCLE_IDs[i], ]
}

Fig2.panSEA <- panSEA::panSEA(GSE12790, all.GSE12790.names, GSEA = FALSE,
                              GSEA.rank.var = rep("Input", length(GSE12790)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              #gmt.features = gmt.features,
                              expression = expression)

# Loading PRISM drug sensitivity AUC scores
# Loading PRISM drug mechanism of action annotations
# Read 85 records
# Read 3655 items
# 123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384851
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# 10
# 11
# 12
# 13
# 14
# 15
# 16
# 17
# 18
# 19
# 20
# 21
# 22
# 23
# 24
# 25
# 26
# 27
# 28
# 29
# 30
# 31
# 32
# 33
# 34
# 35
# 36
# 37
# 38
# 39
# 40
# 41
# 42
# 43
# 44
# 45
# 46
# 47
# 48
# 49
# 50
# 51
# 52
# 53
# 54
# 55
# 56
# 57
# 58
# 59
# 60
# 61
# 62
# 63
# 64
# 65
# 66
# 67
# 68
# 69
# 70
# 71
# 72
# 73
# 74
# 75
# 76
# 77
# 78
# 79
# 80
# 81
# 82
# 83
# 84
# Running DMEA using HCC1806 data
# Running correlations and regressions...
# Error in `[.data.frame`(all.data.corr, , 3:ncol(all.data.corr)) : 
#   undefined columns selected

# save results
fig.folder <- "Fig2"
savePanSEA(Fig2.panSEA, fig.folder)
#Fig2.panSEA <- readRDS(paste0("panSEA_", fig.folder, ".rds"))
Fig2.panSEA <- NULL # make space to process next analysis

#### Step 5b: GSE31625 ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.GSE31625.names)) {
  gmt.features[[all.GSE31625.names[i]]] <- gene.gmt
  expression[[all.GSE31625.names[i]]] <- RNA.df[
    RNA.df$CCLE_ID != all.GSE31625.CCLE_IDs[i], ]
}

Fig3.panSEA <- panSEA::panSEA(GSE31625, all.GSE31625.names, 
                              GSEA.rank.var = rep("Input", length(GSE31625)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)
# save results
fig.folder <- "Fig3"
savePanSEA(Fig3.panSEA, fig.folder)
#Fig3.panSEA <- readRDS(paste0("panSEA_", fig.folder, ".rds"))
Fig3.panSEA <- NULL # make space to process next analysis

#### Step 5c: Coldren ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.Coldren.names)) {
  gmt.features[[all.Coldren.names[i]]] <- gene.gmt
  expression[[all.Coldren.names[i]]] <- RNA.df[
    RNA.df$CCLE_ID != all.Coldren.CCLE_IDs[i], ]
}

Fig4.panSEA <- panSEA::panSEA(Coldren, all.Coldren.names, 
                              GSEA.rank.var = rep("Input", length(Coldren)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)
# save results
fig.folder <- "Fig4"
savePanSEA(Fig4.panSEA, fig.folder)
#Fig4.panSEA <- readRDS(paste0("panSEA_", fig.folder, ".rds"))
Fig4.panSEA <- NULL # make space to process next analysis

#### Step 5d: Coldren proteomics ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.Coldren.prot.names)) {
  gmt.features[[all.Coldren.prot.names[i]]] <- gene.gmt
  expression[[all.Coldren.prot.names[i]]] <- prot.df[
    prot.df$CCLE_ID != all.Coldren.prot.CCLE_IDs[i], ]
}

Fig5.panSEA <- panSEA::panSEA(Coldren.prot.list, all.Coldren.prot.names, 
                              GSEA.rank.var = rep("Input", 
                                                  length(Coldren.prot.list)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)
# save results
fig.folder <- "Fig5"
savePanSEA(Fig5.panSEA, fig.folder)
#Fig5.panSEA <- readRDS(paste0("panSEA_", fig.folder, ".rds"))
Fig5.panSEA <- NULL # make space to process next analysis

#### Step 5e: Coldren metabolomics ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.Coldren.met.names)) {
  gmt.features[[all.Coldren.met.names[i]]] <- gene.gmt
  expression[[all.Coldren.met.names[i]]] <- met.df[
    met.df$CCLE_ID != all.Coldren.met.CCLE_IDs[i], ]
}

Fig6.panSEA <- panSEA::panSEA(Coldren.met.list, all.Coldren.met.names, 
                              GSEA.rank.var = rep("Input", 
                                                  length(Coldren.met.list)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)
# save results
fig.folder <- "Fig6"
savePanSEA(Fig6.panSEA, fig.folder)
#Fig6.panSEA <- readRDS(paste0("panSEA_", fig.folder, ".rds"))
Fig6.panSEA <- NULL # make space to process next analysis

#### Step 5f: Coldren but only querying non-small-cell lung cancer cell line data ####
NSCLC.samples <- cell.line.info[cell.line.info$lineage_subtype == "NSCLC", ] # 273 CCLE
NSCLC.RNA.df <- RNA.df[RNA.df$CCLE_ID %in% NSCLC.samples$CCLE_Name, ] # 63 CCLE
NSCLC.prot.df <- prot.df[prot.df$CCLE_ID %in% NSCLC.samples$CCLE_Name, ] # 78 CCLE
NSCLC.met.df <- met.df[met.df$CCLE_ID %in% NSCLC.samples$CCLE_Name, ] # 182 CCLE
NSCLC.prot.df.noNA <- NSCLC.prot.df[, colSums(is.na(NSCLC.prot.df)) == 0] # 5565 gene names

gmt.features <- list()
expression <- list()
for (i in 1:length(all.Coldren.names)) {
  gmt.features[[all.Coldren.names[i]]] <- gene.gmt
  expression[[all.Coldren.names[i]]] <- NSCLC.RNA.df[
    NSCLC.RNA.df$CCLE_ID != all.Coldren.CCLE_IDs[i], ]
}

Fig7.panSEA <- panSEA::panSEA(Coldren, all.Coldren.names, 
                              GSEA.rank.var = rep("Input", length(Coldren)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)

# save results
setwd(base.path)
savePanSEA(Fig7.panSEA, "Fig5")
Fig7.panSEA <- NULL # make space to process next analysis

#### Step 5g: Coldren proteomics but only querying non-small-cell lung cancer cell line data ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.Coldren.prot.names)) {
  gmt.features[[all.Coldren.prot.names[i]]] <- gene.gmt
  expression[[all.Coldren.prot.names[i]]] <- NSCLC.prot.df.noNA[
    NSCLC.prot.df.noNA$CCLE_ID != all.Coldren.prot.CCLE_IDs[i], ]
}

Fig8.panSEA <- panSEA::panSEA(Coldren.prot.list, all.Coldren.prot.names, 
                              GSEA.rank.var = rep("Input", 
                                                  length(Coldren.prot.list)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)

# save results
setwd(base.path)
savePanSEA(Fig8.panSEA, "Fig8")
Fig8.panSEA <- NULL # make space to process next analysis

#### Step 5h: Coldren proteomics but only querying non-small-cell lung cancer cell line data ####
gmt.features <- list()
expression <- list()
for (i in 1:length(all.Coldren.met.names)) {
  gmt.features[[all.Coldren.met.names[i]]] <- gene.gmt
  expression[[all.Coldren.met.names[i]]] <- NSCLC.met.df[
    NSCLC.met.df$CCLE_ID != all.Coldren.met.CCLE_IDs[i], ]
}

Fig8.panSEA <- panSEA::panSEA(Coldren.met.list, all.Coldren.met.names, 
                              GSEA.rank.var = rep("Input", 
                                                  length(Coldren.met.list)),
                              DMEA.type = "cell_corr", 
                              group.names = group.names,
                              group.samples = group.samples, 
                              gmt.features = gmt.features,
                              expression = expression)

# save results
setwd(base.path)
savePanSEA(Fig8.panSEA, "Fig8")
Fig8.panSEA <- NULL # make space to process next analysis
