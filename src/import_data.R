

if(!exists("chip_index")){
  chip_index <- read.csv("~/projects/supergenes/results/tables/probes_chip_index", header=TRUE, row.names=1)
}



if (!exists("DMR_table")) {
  require(readr)
  print("Loading LUSC DMR metadata...")
  DMR_table <- as.data.frame(read_delim("~/projects/supergenes/data/dataAnais/TCGA-LUSC_cancer_DMRs.tsv", ",", escape_double = FALSE, trim_ws = TRUE))
  
}

if (!exists("cnv_lusc")) {
  print("Loading LUSC copy number variations (cnv) data...")
  cnv_lusc = readRDS("~/projects/tcga_studies/study_TCGA-LUSC_cnv.rds")
}

if (!exists("platform")){
  require(readr)
  platform <- read.csv("~/projects/supergenes/results/tables/biological_regions.csv", row.names=1)
}

if (!exists("cgi_pf")){
  require(readr)
  cgi_pf <- read.csv("~/projects/supergenes/results/tables/cgi_pf.csv", row.names=1)
}


if (!exists("cgi_coordinates")){
  require(readr)
  cgi_coordinates <- read.csv("~/projects/supergenes/results/tables/CGI_coordinates.csv", row.names=1)
}



if (!exists("trscr_lusc")) {
  print("Loading LUSC transcription data...")
  trscr_lusc = readRDS("~/projects/tcga_studies/study_TCGA-LUSC_trscr.rds")
}

if (!exists("meth_lusc")) {
  print("Loading LUSC methylation data...")
  meth_lusc = readRDS("~/projects/tcga_studies/study_TCGA-LUSC_meth.rds")
  meth_tumoral <- meth_lusc$data[,rownames(meth_lusc$exp_grp[which(meth_lusc$exp_grp[,"tissue_status"]=="tumoral"),])]
  meth_normal <- meth_lusc$data[,rownames(meth_lusc$exp_grp[which(meth_lusc$exp_grp[,"tissue_status"]=="normal"),])]
}
if (!exists("meth_tumoral")){
  print("Spliting methylation data (Tumoral, healthy, differential...")
  meth_tumoral <- meth_lusc$data[,rownames(meth_lusc$exp_grp[which(meth_lusc$exp_grp[,"tissue_status"]=="tumoral"),])]
  meth_normal <- meth_lusc$data[,rownames(meth_lusc$exp_grp[which(meth_lusc$exp_grp[,"tissue_status"]=="normal"),])]
  mean_values_normal_tissues <- apply(meth_normal,1,mean,na.rm=T)
  meth_diff <- t(sapply(1:nrow(meth_tumoral), function(i) meth_tumoral[i,] - mean_values_normal_tissues[i]))
  rownames(meth_diff)<-rownames(meth_lusc$data)
}

if(!exists("penda_superup_deregulated")) {
  print("Loading superexpressed genes list...")
  penda_superup_deregulated = as.data.frame(readxl::read_excel("~/projects/supergenes/data/tables_penda.xlsx"))
  rownames(penda_superup_deregulated) <- penda_superup_deregulated[,1]; penda_superup_deregulated <- penda_superup_deregulated[,-1]
  
  for (i in 1:dim(penda_superup_deregulated)[1]){
    if (is.element(i,c(1:3,178:dim(penda_superup_deregulated)[1])))
    {penda_superup_deregulated[i,1] <- 1}
    else {penda_superup_deregulated[i,1] <- 0}
    
    if (i <= 177)
    {penda_superup_deregulated[i,2] <- 1}
    else {penda_superup_deregulated[i,2] <- 0}
  }
  penda_superup_deregulated <- rownames(penda_superup_deregulated[penda_superup_deregulated[,2] == 1, ])
}  
if(!exists("penda_superdown_deregulated")){
  print("Loading underexpressed genes list...")
  penda_superdown_deregulated = as.data.frame(readxl::read_excel("~/projects/supergenes/data/tables_penda.xlsx",sheet = 2))
  rownames(penda_superdown_deregulated) <- penda_superdown_deregulated[,1];
  penda_superdown_deregulated <- penda_superdown_deregulated[,-1]
  
  for (i in 1:dim(penda_superdown_deregulated)[1]){
    if (is.element(i,c(1:123,385:dim(penda_superdown_deregulated)[1]))) 
    {penda_superdown_deregulated[i,1] <- 1}
    else {penda_superdown_deregulated[i,1]<- 0}
    
    if(i <= 384)
    {penda_superdown_deregulated[i,2] <- 1}
    else {penda_superdown_deregulated[i,2]<- 0}
    
  }
  
  penda_superdown_deregulated <- rownames(penda_superdown_deregulated[penda_superdown_deregulated[,2] == 1, ])
}


if(!exists("penda_superconserved")){
  print("Loading conserved genes list...")
  penda_superconserved = as.data.frame(readxl::read_excel("~/projects/supergenes/data/tables_penda.xlsx",sheet = 3))
  rownames(penda_superconserved) <- penda_superconserved[,1]; penda_superconserved <- penda_superconserved[,-1]
  
  for (i in 1:dim(penda_superconserved)[1]){
    if (is.element(i,c(1:33,80:dim(penda_superconserved)[1])))
    {penda_superconserved[i,1] <- 1}
    else {penda_superconserved[i,1] <- 0}
    
    if(i <= 79)
    {penda_superconserved[i,2] <- 1}
    else {penda_superconserved[i,2] <- 0}
  }
  penda_superconserved <- rownames(penda_superconserved[penda_superconserved[,2] == 1, ])
  
}



if (!exists("epic")) {
  print("Loading Epic metadata .")
  epic_orig = readRDS("~/projects/datashare/platforms/EPIC.hg38.manifest.gencode.v22.rds")
  # epic_orig = readRDS("~/projects/datashare/platforms/EPIC.hg38.manifest.rds")
  epic = as.data.frame(epic_orig)
  rownames(epic) = names(epic_orig)
  head(epic)
  pf_chr_colname = "seqnames"
  pf_pos_colname = "start"
  epic = epic[order(epic[[pf_chr_colname]], epic[[pf_pos_colname]]),]
  rm(epic_orig)
}


meth_tumoral <- meth_lusc$data[,rownames(meth_lusc$exp_grp[which(meth_lusc$exp_grp[,"tissue_status"]=="tumoral"),])]
meth_normal <- meth_lusc$data[,rownames(meth_lusc$exp_grp[which(meth_lusc$exp_grp[,"tissue_status"]=="normal"),])]

source("~/projects/supergenes/src/fun_preprocess.R")
source("~/projects/supergenes/src/fun_plots.R")
source("~/projects/supergenes/src/fun_calculus.R")
gc()