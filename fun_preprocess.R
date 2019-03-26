###################################  Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################


get_features <- function(targeted_genes, index_pathology = 2, study, up_str = 2500, dwn_str = 2500, nb_probe_min = 1){   #### In my uses-case, 1 for LUAD & 2 for LUSC 
  print("Creating a list of features...")
  targeted_genes = targeted_genes[targeted_genes[,index_pathology]==1,]
  features = study$platform[intersect(rownames(study$platform),rownames(targeted_genes)),]
  
  TSS = c()
  for (i in 1:dim(features)[1]){
    if (features$strand[i] == "-"){TSS[i] <- features$tx_end[i]}
    else{TSS[i] <- features$tx_start[i]}
  }
  
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
  feat_indexed_probes = epimedtools::monitored_apply(features, 1, function(gene) {
    # gene = features[1,]
    # print(gene)
    chr = gene[[1]]
    meth_platform = chrs_indexed_epic[[chr]]
    ret = dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str) 
    return(ret)
  })
  
  
  features$nb_epic_probes = sapply(feat_indexed_probes[rownames(features)], length)
  features = cbind(features,TSS)
  sub_features =  features[features$nb_epic_probes >= nb_probe_min,]
  
  warning(paste("\n",nrow(features)-nrow(sub_features),"genes were removed due to the followings selection parameters :","\n",
                "window downstream the TSS :",dwn_str,"\n",
                "window upstream the TSS : ", up_str, "\n",
                "Min number of probes :", nb_probe_min, sep = " "))
  
  ret = list("features" = sub_features, "probes" = feat_indexed_probes)
  return(ret)
  
}

#Genes selected catalog function ( no calculus only visualisation)


catalog =  function(selected_features = features){
  
  available_genes = list(catalog="catalog", n_genes = "n_genes")
  
  available_genes[["catalog"]] = noquote(rownames(selected_features$features))
  available_genes [["n_genes"]] = noquote(paste("Genes selected :", length(available_genes[[1]]),sep = " "))
  warning("Note that quotes were removed for clarity, please select your genes using it for the plotting function.")
  return(available_genes)
  
}



###########################  Creation of bins ######################

get_bins = function(selected_gene = "TG", binwidth = bw, nbins =20){
  window = c(features$features[selected_gene,"TSS"]-(bw*nbins),features$features[selected_gene,"TSS"]+(bw*nbins))
  window_width = window[2]-window[1]
  bins_coordinates = seq(window[1],window[2],bw)
  tmp = meth_lusc$platform[selected_probes,]
  
  
  
  
  
}












