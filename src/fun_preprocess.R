###################################  Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################


get_features <- function(targeted_genes, study, up_str = 2500, dwn_str = 2500, nb_probe_min = 1) { 
  print("Creating a list of features...")
  features <- study$platform[intersect(rownames(study$platform), targeted_genes), ]

  TSS <- c()
  for (i in 1:dim(features)[1]) {
    if (features$strand[i] == "-") {
      TSS[i] <- features$tx_end[i]
    }
    else {
      TSS[i] <- features$tx_start[i]
    }
  }

  print("Indexing probe by features")
  # params
  ## index meth probes by chr
  pf_chr_colname <- "seqnames"
  pf_pos_colname <- "start"
  chrs <- unique(features[, 1])
  chrs_indexed_epic <- lapply(chrs, function(chr) {
    print(chr)
    idx <- rownames(epic)[epic[[pf_chr_colname]] %in% chr]
    ret <- epic[idx, ]
    return(ret)
  })
  names(chrs_indexed_epic) <- chrs

  ## index probes by gene name
  feat_indexed_probes <<- epimedtools::monitored_apply(features, 1, function(gene) {
    # gene = features[1,]
    # print(gene)
    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[chr]]
    tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]
    # here compute bin 1 to 6 using sub_epic
    # do_biens
    # ret = list(sub_epic=sub_epic, bin1=..., bin2=..., bin3=..., bin4=..., bin5=..., bin6=...)
    return(tmp_probes)
  })


  features$nb_epic_probes <- sapply(feat_indexed_probes[rownames(features)], length)
  features <- cbind(features, TSS)
  sub_features <- features[features$nb_epic_probes >= nb_probe_min, ]

  warning(paste("\n", nrow(features) - nrow(sub_features), "genes were removed due to the followings selection parameters :", "\n",
    "window downstream the TSS :", dwn_str, "\n",
    "window upstream the TSS : ", up_str, "\n",
    "Min number of probes :", nb_probe_min,
    sep = " "
  ))

  return(sub_features)
}




############################### get_features_bins : same as prior + indexing probes per bins #####################

get_features_bins <- function(targeted_genes, study = trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1, ...) { #### In my uses-case, 1 for LUAD & 2 for LUSC
  print("Creating a list of features...")
  features <- study$platform[intersect(rownames(study$platform), targeted_genes), ]
  
  
  
  
  print("Indexing probes by features and by bins")
  
  pf_chr_colname <- "seqnames"
  pf_pos_colname <- "start"
  chrs <- unique(features[, 1])
  chrs_indexed_epic <- lapply(chrs, function(chr) {
    print(chr)
    idx <- rownames(epic)[epic[[pf_chr_colname]] %in% chr]
    ret <- epic[idx, ]
    return(ret)
  })
  names(chrs_indexed_epic) <- chrs
  
  
  
  
  feat_indexed_probes_bin <- epimedtools::monitored_apply(features, 1, function(gene) {
    
    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[chr]]
    tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) - 2500
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) + 1000
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 1500
    bin1 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) - 1000
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) + 500
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin2 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) - 500
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) + 0
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin3 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) + 0
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) - 500
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin4 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) + 500
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) - 1000
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin5 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) + 1000
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) - 2500
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 1500
    bin6 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    whole_feature <- tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    
    
    ret <- list(sub_epic = sub_epic, whole_feature = whole_feature, bin1 = bin1, bin2 = bin2, bin3 = bin3, bin4 = bin4, bin5 = bin5, bin6 = bin6)
    return(ret)
    
    
  })
  
  
  
  return(feat_indexed_probes_bin)
}

############catalog##################

catalog <- function(selected_features = features) {
  available_genes <- list(catalog = "catalog", n_genes = "n_genes")
  
  available_genes[["catalog"]] <- noquote(rownames(selected_features))
  available_genes [["n_genes"]] <- noquote(paste("Genes selected :", length(available_genes[[1]]), sep = " "))
  warning("Note that quotes were removed for clarity, please select your genes using it for the plotting function.")
  return(available_genes)
}


############## reduce_rows #################

reduce_rows <- function(tmp_meth_data, map, indicator_func2 = mean, ...) {
  dim(tmp_meth_data)
  meth_by_tissues_by_feature <- epimedtools::monitored_apply(mod = 10, t(t(names(map))), 1, function(f) {
    # print(f)
    # f = features[1,]
    # f = features["HACD4",]
    probe_idx <- intersect(rownames(tmp_meth_data), map[[f]])
    if (length(probe_idx) == 0) {
      # print(f)
      tmp_meth_by_tissues_by_feature <- rep(NA, ncol(tmp_meth_data))
      names(tmp_meth_by_tissues_by_feature) <- colnames(tmp_meth_data)
    } else if (length(probe_idx) > 1) {
      tmp_meth_by_tissues_by_feature <- apply(tmp_meth_data[probe_idx, ], 2, indicator_func2, ...)
    } else {
      tmp_meth_by_tissues_by_feature <- sapply(tmp_meth_data[probe_idx, ], indicator_func2, ...)
    }
    return(tmp_meth_by_tissues_by_feature)
  })
  colnames(meth_by_tissues_by_feature) <- names(map)
  meth_by_tissues_by_feature <- t(meth_by_tissues_by_feature)
  return(meth_by_tissues_by_feature)
}





############### get a one level map for genes per bin with a form such as 


get_binmap <- function(feat_indexed_probes_bin = feat_indexed_probes_bin, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"))  {
  bin_indexed_probes = lapply(feat_indexed_probes_bin, function(g){
    g[binlist]
  })
  bin_indexed_probes = unlist(bin_indexed_probes, recursive = FALSE)
}

####################### get means per bins per genes #############################

#means_per_bins_per_genes_per_patient = reduce_rows(meth_lusc$data,binmap,mean,na.rm=T)

subset_vals_per_bins <- function(values_per_patient = means_per_bins_per_genes_per_patient, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"), fun = mean, ...){
  values_per_bins <- apply(values_per_patient,1,fun, na.rm=T) ### get mean per lines (gene_bin)
  ret <- list()
  
  vals <- sapply(binlist, function(bin){
    ret[[bin]] <- values_per_bins[endsWith(names(values_per_bins),bin)]    
    return(ret)
    
    }
  )
  
  
  tmp = reduce_rows(meth_lusc$data, feat_indexed_probes, mean, na.rm=T)
  overall = apply(tmp,1,fun,na.rm=T)
  
  
  vals_per_genes <- do.call(cbind,vals)
  vals_per_genes <- cbind(vals_per_genes,overall)
  colnames(vals_per_genes) <-c(binlist,"overall")
  rownames(vals_per_genes) <- names(feat_indexed_probes_bin)
  
  return(vals_per_genes)
  
}

#################heatmap######################

meth_heatmap <- function(data = means, dendrogram = "none", Rowv = NULL, Colv = NULL, ...){
  require(gplots)
  
  data = data[order(data[,7]),]
  data = data[,-7]
  
  
  
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(20)
  
  foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, trace="none", col=cols, main=paste0("Mean of values (", nrow(data), " genes x ", ncol(data), "bins"), mar=c(10,5), useRaster=TRUE)
  
  return(foo)
}
