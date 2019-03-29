get_features <- function(targeted_genes, index_pathology = 2, study = trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1){   #### In my uses-case, 1 for LUAD & 2 for LUSC 
print("Creating a list of features...")
targeted_genes = targeted_genes[targeted_genes[,index_pathology]==1,]
features = study$platform[intersect(rownames(study$platform),rownames(targeted_genes)),]



  
  print("Indexing probe by features")
  # params
  ## index meth probes by chr
  pf_chr_colname = "seqnames"
  pf_pos_colname = "start"
  chrs = unique(features[,1])
  chrs_indexed_epic = lapply(chrs, function(chr) {
    print(chr)
    idx = rownames(epic)[epic[[pf_chr_colname]] %in% chr]  
    ret = epic[idx,]
    return(ret)
  })
  names(chrs_indexed_epic) = chrs
  
 
  
  ## index probes by gene name
  feat_indexed_probes_bin <<-epimedtools::monitored_apply(features, 1, function(gene){
    # gene = features[1,]
    # print(gene)
    chr = gene[[1]]
    meth_platform = chrs_indexed_epic[[chr]]
    tmp_probes = dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    sub_epic = meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]
    # here compute bin 1 to 6 using sub_epic
    # do_bins()
    tmp_gene = gene
    tmp_gene[[6]] = "+"
    if (gene[[6]]=="+") {
      tmp_gene[[2]] = as.numeric(tmp_gene[[2]]) - 2500
    } else {
      tmp_gene[[2]] = as.numeric(tmp_gene[[3]]) + 1000
    }
    tmp_u_bin = 0
    tmp_d_bin = 1500
    bin1 = dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
  
    
    tmp_gene = gene
    tmp_gene[[6]] = "+"
    if (gene[[6]] =="+") {
      tmp_gene[[2]] = as.numeric(tmp_gene[[2]]) - 1000
    } else {
      tmp_gene[[2]] = as.numeric(tmp_gene[[3]]) + 500
    }
    tmp_u_bin = 0
    tmp_d_bin = 500
    bin2 = dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene = gene
    tmp_gene[[6]] = "+"
    if (gene[[6]]=="+") {
      tmp_gene[[2]] = as.numeric(tmp_gene[[2]]) - 500
    } else {
      tmp_gene[[2]] = as.numeric(tmp_gene[[3]]) + 0
    }
    tmp_u_bin = 0
    tmp_d_bin = 500
    bin3 = dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene = gene
    tmp_gene[[6]] = "+"
    if (gene[[6]]=="+") {
      tmp_gene[[2]] = as.numeric(tmp_gene[[2]]) + 0
    } else {
      tmp_gene[[2]] = as.numeric(tmp_gene[[3]]) - 500
    }
    tmp_u_bin = 0
    tmp_d_bin = 500
    bin4 = dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene = gene
    tmp_gene[[6]] = "+"
    if (gene[[6]]=="+") {
      tmp_gene[[2]] = as.numeric(tmp_gene[[2]]) + 500
    } else {
      tmp_gene[[2]] = as.numeric(tmp_gene[[3]]) - 1000
    }
    tmp_u_bin = 0
    tmp_d_bin = 500
    bin5 = dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene = gene
    tmp_gene[[6]] = "+"
    if (gene[[6]]=="+") {
      tmp_gene[[2]] = as.numeric(tmp_gene[[2]]) +1000
    } else {
      tmp_gene[[2]] = as.numeric(tmp_gene[[3]]) - 2500
    }
    tmp_u_bin = 0
    tmp_d_bin = 1500
    bin6 = dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    ret = list(sub_epic=sub_epic, bin1=bin1, bin2=bin2, bin3=bin3, bin4=bin4, bin5=bin5, bin6=bin6)
    return(ret)
    
  })
  
  
  features$nb_epic_probes = sapply(feat_indexed_probes_bin[rownames(features)], length)
  sub_features =  features[features$nb_epic_probes >= nb_probe_min,]
  
  print(paste(nrow(features)-nrow(sub_features),"genes were removed due to the followings selection parameters :",
                "window downstream the TSS :",dwn_str,
                "window upstream the TSS : ", up_str,
                "Min number of probes :", nb_probe_min, sep = " "))
  
  
  return(sub_features)
}
