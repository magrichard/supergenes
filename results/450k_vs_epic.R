############### sondes epic / données meth############

if (!exists("meth_lusc")) {
  print("Loading LUSC methylation data...")
  meth_lusc = readRDS("~/projects/supergenes/data/tcga_studies/study_TCGA-LUSC_meth.rds")
}

if (!exists("epic")) {
  print("Loading Epic metadata .")
  epic_orig = readRDS("~/projects/supergenes/data/EPIC.hg38.manifest.gencode.v22.rds")
  # epic_orig = readRDS("~/projects/datashare/platforms/EPIC.hg38.manifest.rds")
  epic = as.data.frame(epic_orig)
  rownames(epic) = names(epic_orig)
  head(epic)
  pf_chr_colname = "seqnames"
  pf_pos_colname = "start"
  epic = epic[order(epic[[pf_chr_colname]], epic[[pf_pos_colname]]),]
}


meth_data_2 <- subset(meth_lusc$data, !(rownames(meth_lusc$data) %in% rownames(epic)))
chr_2 <- meth_lusc$platform[rownames(meth_data_2),1]
length(chr_2)
table(chr_2)
prop.table(table(chr_2))
barplot(table(chr_2))