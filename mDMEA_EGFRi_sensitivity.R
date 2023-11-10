# Use BeatAML data for DMEA
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
# Created: 2023-10-10
# Last modified: 2023-10-17

#### Step 1. Prepare environment ####
## define custom functions
# set illegal filename characters
illegal.chars <- c("#","<","%",">","!","`","&","'","=","}","/",":","@") # source: https://www.mtu.edu/umc/services/websites/writing/characters-avoid/; would be nice to protect against " and \ too
illegal.chars.need.brackets <- c("$","+","*","|","{","?") # for these chars, gsub needs brackets but str_contains can't have brackets

as.filename <- function(moa.name, 
                        illegal.chars = c("#","<","%",">","!","`","&","'","=","}","/",":","@"), 
                        illegal.chars.need.brackets = c("$","+","*","|","{","?")){
  # if moa.name contains illegal file name characters, replace these characters with "-" or "_"
  if(sjmisc::str_contains(moa.name, illegal.chars, logic="or")){
    for(j in 1:length(illegal.chars)){
      moa.name <- gsub(illegal.chars[j], "-", moa.name)
    }
  }
  if(sjmisc::str_contains(moa.name, illegal.chars.need.brackets, logic="or")){
    for(j in 1:length(illegal.chars.need.brackets)){
      moa.name <- gsub(paste0("[",illegal.chars.need.brackets[j],"]"), "-", moa.name)
    }
  }
  if(sjmisc::str_contains(moa.name, " ")){
    moa.name <- gsub(" ", "_", moa.name)
  }
  return(moa.name)
}

## install packages from GitHub
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

# DMEA
if (!require("DMEA", quietly = TRUE))
  devtools::install_github("BelindaBGarana/DMEA")

## install packages non-default repos
if (!require("synapser", quietly = TRUE))
  install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

## install CRAN packages
# reshape2
if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2")

# ggplot2
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

## load packages into environment
library(DMEA);library(synapser);library(utils);library(reshape2);library(ggplot2);


#### Step 2. Import & format BeatAML data ####
## log-in to Synapse
#synapser::synLogin() # use your Synapse login credentials after disconnecting from vpn

## drug MOA annotations ## REPLACE with GitHub link when repo is public
moa.df <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

## differential expression data from PTRC2 global proteomics of AML cell lines (project synID: syn39679334)
exps <- c("12_&_18", "18")
diff.syn <- list(synapser::synGet('syn51198695'), # exps 12 & 18 combined
                 synapser::synGet('syn51514584')) # exp 18
# exp 12 global data is syn22156807 but not analyzed for diffexp yet; exp 12 phospho data is syn22156809
# exp 18 phospho diffexp is syn51514583
diff.df <- list()
for(i in 1:length(exps)){
  exp.name <- paste0("exp_", exps[i])
  if (exps[i] == "18") {
    diff.df[[exp.name]] <- readxl::read_excel(diff.syn[[i]]$path)
    
    # filter diffexp for all adjusted p-values < 0.05
    diff.df[[exp.name]] <- diff.df[[exp.name]][diff.df[[exp.name]]$limma_adj < 0.05 & 
                                                 diff.df[[exp.name]]$t_test_adj < 0.05 & 
                                                 diff.df[[exp.name]]$welch_adj < 0.05, ]
  } else {
    diff.df[[exp.name]] <- utils::read.table(diff.syn[[i]]$path, sep = "\t", header = TRUE)
    
    # filter diffexp for all adjusted p-values < 0.05
    diff.df[[exp.name]] <- diff.df[[exp.name]][diff.df[[exp.name]]$limma_adj < 0.05 & 
                                                 diff.df[[exp.name]]$t_test_adj < 0.05 & 
                                                 diff.df[[exp.name]]$welch_adj < 0.05, ]
  }
}

## Proteomic data from Synapse (project synID: syn22128879)
# metadata for patients in BeatAML
meta.syn <- synapser::synGet('syn25807733')
meta.df <- utils::read.table(meta.syn$path, sep = "\t", header = TRUE)

# global proteomics from patients in BeatAML
global.syn <- synapser::synGet('syn25714248')
global.df <- utils::read.table(global.syn$path, sep = "\t", header = TRUE)

# drug AUC data from patients in BeatAML
drug.syn <- synapser::synGet('syn51674470')
drug.df <- utils::read.csv(drug.syn$path)

## transcriptomic data from BeatAML waves 1-4: https://biodev.github.io/BeatAML2/
options(timeout = 120) # loading RNAseq data takes >60s (default timeout)
RNAseq.df <- utils::read.table(
  paste0("https://github.com/biodev/beataml2.0_data/raw/main/",
         "beataml_waves1to4_norm_exp_dbgap.txt"), 
  sep = "\t", header = TRUE) # gene names are in column display_label (#2), sample names are colnames 5-end

drug14.df <- utils::read.table(
  paste0("https://github.com/biodev/beataml2.0_data/raw/main/",
         "beataml_probit_curve_fits_v4_dbgap.txt"), 
  sep = "\t", header = TRUE) # sample names for RNAseq.df are in column dbgap_rnaseq_sample (#3), drug names in column Iihibitor (#4), auc values in column auc (#23)

# only keep drug data for samples with RNAseq data
drug14.df <- drug14.df[drug14.df$dbgap_rnaseq_sample != "", ]

## metabolomic data from BeatAML PTRC2 (project synID: syn39679334)
met.syn <- synapser::synGet('syn52224584')
met.df <- readxl::read_excel(met.syn$path)

### format data for DMEA
## Proteomic data from Synapse (project synID: syn22128879)
# format drug.df wide (samples in first column, drug names for rest of columns)
drug.df <- reshape2::dcast(drug.df, sample_id ~ inhibitor, value.var = "auc", fill = NA)

# remove drugs without moa annotations and drug combos
valid.drugs <- names(drug.df)[names(drug.df) %in% moa.df[!is.na(moa.df),]$Drug] # 167 drugs
drug.df <- drug.df[ , c("sample_id", valid.drugs)] # 167 drugs
moa.df <- moa.df[moa.df$Drug %in% names(drug.df)[2:ncol(drug.df)], ]

# change global.df column names from SampleID.abbrev to Barcode.ID to match drug.df
global.ids <- names(global.df)

# remove X and any 0's from start of each column name
# replace SampleID.abbrev with Barcode.ID to match drug.df
for(i in seq_len(length(global.ids))){
  global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  
  if(substring(global.ids[i], 1, 1) == 0){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  }
  
  if(global.ids[i] %in% meta.df$SampleID.abbrev){
    global.ids[i] <- meta.df[meta.df$SampleID.abbrev == global.ids[i], ]$Barcode.ID
  }
}

# replace global.df column names 
names(global.df) <- global.ids

# transpose global.df so that first column is Barcode.ID and rest of columns are gene symbols
global.df <- as.data.frame(t(global.df))

# make first column Barcode.ID
global.df$Barcode.ID <- rownames(global.df)
global.df <- global.df[ , c("Barcode.ID", names(global.df[ , 1:(ncol(global.df)-1)]))]

# make sure columns with sample IDs have same names between global.df and drug.df
names(drug.df)[1] <- names(global.df)[1]

## transcriptomic data from BeatAML waves 1-4: https://biodev.github.io/BeatAML2/
# format drug.df wide (samples in first column, drug names for rest of columns)
drug14.df <- reshape2::dcast(drug14.df, dbgap_rnaseq_sample ~ inhibitor, value.var = "auc", fill = NA)

# remove drugs without moa annotations
valid14.drugs <- names(drug14.df)[names(drug14.df) %in% moa.df[!is.na(moa.df),]$Drug] # 165 drugs ("Arsenic Trioxide" and "F_S_V" from valid.drugs are not in valid14.drugs, 
# but F_S_V does not have moa annotation and Arsenic Trioxide's moa antineoplastic agent isn't used in DMEA)
drug14.df <- drug14.df[ , c("dbgap_rnaseq_sample", valid14.drugs)] # 165 drugs
#moa14.df <- moa.df[moa.df$Drug %in% names(drug14.df)[2:ncol(drug14.df)], ]

# transpose RNAseq.df so that first column is dbgap_rnaseq_sample and rest of columns are gene symbols
RNAseq.df <- RNAseq.df[ , c(2, 5:ncol(RNAseq.df))]
rownames(RNAseq.df) <- RNAseq.df$display_label
RNAseq.df$display_label <- NULL
RNAseq.df <- as.data.frame(t(RNAseq.df))
RNAseq.df$dbgap_rnaseq_sample <- rownames(RNAseq.df)
RNAseq.df <- RNAseq.df[ , c("dbgap_rnaseq_sample", names(RNAseq.df)[1:(ncol(RNAseq.df)-1)])]

#### Step 3. Import & format CCLE data ####

#### Step 5. Run DMEA to identify selectively toxic drug MOAs for each exp, contrast, data set ####
DMEAwd <- "~/OneDrive - PNNL/Documents/PTRC2/DMEA/Using_BeatAML/"
contrasts <- na.omit(unique(c(diff.df[[exps[1]]]$contrast, 
                              diff.df[[exps[2]]]$contrast)))
types <- c("Transcriptomics", "Proteomics", "Metabolomics")
sources <- c("BeatAML", "CCLE")

RNA <- list(RNAseq.df, CCLE.RNAseq.df)
prot <- list(global.df, CCLE.global.df)
met <- list(met.df, CCLE.met.df)

BeatAML.DMEA <- list()
global.DMEA <- list()
met.DMEA <- list()
for (i in 1:length(contrasts)) {
  dir.create(contrasts[i])
  setwd(contrasts[i])
  
  RNAseq.DMEA[[contrasts[i]]] <- list()
  global.DMEA[[contrasts[i]]] <- list()
  met.DMEA[[contrasts[i]]] <- list()
  
  for (j in 1:length(exps)) {
    # check if contrast i was studied for exp j
    if (contrasts[i] %in% unique(diff.df[[exps[j]]]$contrast)) {
      exp.name <- paste0("exp_", exps[i])
      dir.create(exp.name)
      setwd(exp.name)
      
      # filter diff.df for each contrast & global proteomics
      contrast.df <- diff.df[[exp.name]][diff.df[[exp.name]]$contrast == contrasts[i] & 
                                           diff.df[[exp.name]]$data_type == "Global", ]
      
      for (k in 1:length(sources)) {
        ### run DMEA for each contrast
        expr <- list(RNAseq.df, global.df, met.df)
        sig <- list(RNAseq.contrast, global.contrast, met.contrast)
        contrast.DMEA[[exp.name]] <- mDMEA::mDMEA(drug.sensitivity = drug.df, 
                                                  expression = expr, 
                                                  weights = sig, types = types,
                                                  drug.info = moa.df, sep = ", ")
        
        # save CSV files
        utils::write.csv(contrast.DMEA[[exp.name]]$drug.est.df, 
                         paste0("DMEA_correlations_", Sys.Date(), ".csv"), 
                         row.names = FALSE)
        utils::write.csv(contrast.DMEA[[exp.name]]$results, 
                         paste0("DMEA_results_", Sys.Date(), ".csv"), 
                         row.names = FALSE)
        
        # save plots as PDFs
        # save DMEA results locally
        setwd(DMEAwd)
        for (m in 1:length(types)) {
          dir.create(types[m])
          setwd(types[m]) 
        }
        
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
      }
    }
    
    RNAseq.DMEA[[contrasts[i]]] <- RNAseq.contrast.DMEA
    global.DMEA[[contrasts[i]]] <- global.contrast.DMEA
    met.DMEA[[contrasts[i]]] <- met.contrast.DMEA
  }
}

## save R result objects for any future analyses
for (i in 1:length(sources)) {
  setwd(paste0(DMEAwd), sources[i])
  saveRDS(global.DMEA, paste0("DMEA_", sources[i], "_", Sys.Date(), ".rds"))
}

#### Step 6. Compile DMEA results ####
# across exps

# across contrasts

# across exps & contrasts

# across sources (BeatAML vs. CCLE)

# across contrasts & sources

# across exps & sources

# across exps, sources & contrasts
