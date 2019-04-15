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





################# Index regions per features################
get_features_regions <- function(features_list = features, platform_regions = platform, epimed_database = epic, indexed_probes = feat_indexed_probes ){
  
  regions_name <-unique(as.character(platform[,"genomic_feature"]))
  
  feat_indexed_probes_regions <- lapply(1:nrow(features_list), function(index){
    
    
    gene <- names(feat_indexed_probes[index])
    tmp_probes <- feat_indexed_probes[[gene]]
    epic_tmp <- epic[tmp_probes,c("start","end")]
    sub_epic <- epic_tmp[intersect(rownames(epic_tmp), rownames(meth_lusc$data)),]
    pf_gene <- platform[which(platform[,"gene"]==gene),]
    
    
    pf_regions = lapply(regions_name, function(region){
      regions<-pf_gene[which(pf_gene[,"genomic_feature"]==region),]
      return(regions)
    })
    
    names(pf_regions) <- regions_name
    
    per_regions_type_coordinates <-sapply(pf_regions, function(pf){
      if (nrow(pf)>=1){
        per_regions_coordinates = apply(pf,1,function(region){
          range = c(region["start"],region["end"])
        })
      }
    })
    
    names(per_regions_type_coordinates) <- regions_name
    
    INTER = vector()
    P2000 = vector()
    UTR5 = vector()
    INTRON = vector()
    CDS = vector()
    UTR3 = vector()
    
    for (i in 1:nrow(sub_epic)){
      
      if (!is.null(per_regions_type_coordinates$INTER))
        if(sum(apply(per_regions_type_coordinates$INTER,2, function(window){
          sub_epic[i,"start"] > window["start"] & sub_epic[i,"end"] < window["end"] }),na.rm=T)>=1)
        {INTER[i]<-rownames(sub_epic[i,])}
      
      
      
      if (!is.null(per_regions_type_coordinates$P2000))
        if(sum(apply(per_regions_type_coordinates$P2000,2, function(window){
          sub_epic[i,"start"] > window["start"] & sub_epic[i,"end"] < window["end"] }),na.rm=T)>=1)
        {P2000[i]<-rownames(sub_epic[i,])}
      
      
      
      if (!is.null(per_regions_type_coordinates$UTR5))
        if(sum(apply(per_regions_type_coordinates$UTR5,2, function(window){
          sub_epic[i,"start"] > window["start"] & sub_epic[i,"end"] < window["end"] }),na.rm=T)>=1)
        {UTR5[i]<-rownames(sub_epic[i,])}
      
      
      
      if (!is.null(per_regions_type_coordinates$INTRON))
        if(sum(apply(per_regions_type_coordinates$INTRON,2, function(window){
          sub_epic[i,"start"] > window["start"] & sub_epic[i,"end"] < window["end"] }),na.rm=T)>=1)
        {INTRON[i]<-rownames(sub_epic[i,])}
      
      
      if (!is.null(per_regions_type_coordinates$CDS))
        if(sum(apply(per_regions_type_coordinates$CDS,2, function(window){
          sub_epic[i,"start"] > window["start"] & sub_epic[i,"end"] < window["end"] }),na.rm=T)>=1) 
        {CDS[i]<-rownames(sub_epic[i,])}
      
      if (!is.null(per_regions_type_coordinates$UTR3))
        if(sum(apply(per_regions_type_coordinates$UTR3,2, function(window){
          sub_epic[i,"start"] > window["start"] & sub_epic[i,"end"] < window["end"] }),na.rm=T)>=1)
        {UTR3[i]<-rownames(sub_epic[i,])}
      
    }
    
    ret <- list(sub_epic = sub_epic, INTER = INTER, P2000 = P2000, UTR5 = UTR5, CDS = CDS, INTRON = INTRON, UTR3 = UTR3)
    return(ret)
    
  })
  
  names(feat_indexed_probes_regions) <- rownames(features_list)
  return(feat_indexed_probes_regions)
  
}

################ Index DMR per features##################

get_indexed_DMRs <- function(DMRtable = DMR_table, features_list = features, regions_coordinates = platform, treshold = 7500){
  
  DMRtable<- DMRtable[order(DMRtable[,"DMR_id"]),]
  
  print("Extractig regions coordinates per available genes...")
  
  pf<- lapply(rownames(features_list), function(gene){
    regions_coordinates[which(regions_coordinates[,"gene"]==gene),]
  })
  names(pf)<- rownames(features_list)
  
  
  
  print("Indexing DMRs per genes...")
  
  feat_indexed_DMRs <- lapply(1:length(rownames(features_list)), function(gene){
    print(paste0("Indexing DMRs for ", rownames(features_list)[gene] ,"..."))
    DMR_in_range <- lapply(1:nrow(DMRtable),function(index){
      DMR <- DMRtable[index,]
      
      
      if(as.numeric(DMR[["start"]]) > features_list[gene,"TSS"]){
        foo <- abs(as.numeric(DMR[["start"]])-features_list[gene,"TSS"])
        return(foo)
      }
      
      if(as.numeric(DMR[["start"]]) < features_list[gene,"TSS"]){
        foo <- abs(features_list[gene,"TSS"]-as.numeric(DMR[["end"]]))
        return(foo)
      }
    })
    
    names(DMR_in_range)<-DMRtable[["DMR_id"]]
    DMR_candidates <- unlist(DMR_in_range)[DMR_in_range <= treshold]
    return(DMR_candidates)
  })
  names(feat_indexed_DMRs)<- rownames(features_list)
  
  
  return(feat_indexed_DMRs)
  
  
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


get_binmap <- function(map_to_reduce = feat_indexed_probes_bin, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"))  {
  bin_indexed_probes = lapply(map_to_reduce, function(g){
    g[binlist]
  })
  bin_indexed_probes = unlist(bin_indexed_probes, recursive = FALSE)
}

####################### get means per bins per genes #############################

#means_per_bins_per_genes_per_patient = reduce_rows(meth_lusc$data,binmap,mean,na.rm=T)

subset_vals_per_bins<-function(data = meth_lusc$data,
                               values_per_patient = means_per_regions_per_genes_per_patient,
                               binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"),
                               fun = mean,
                               names = feat_indexed_probes, ...){
  
  values_per_bins <- apply(values_per_patient,1,mean , na.rm=T) 
  ret <- list()
  
  vals <- sapply(binlist, function(bin){
    ret[[bin]] <- values_per_bins[endsWith(names(values_per_bins),bin)]    
    return(ret)
    
  }
  )
  
  
  tmp = reduce_rows(data, feat_indexed_probes, fun, na.rm=T) ### note that you need feat_indexed_probes loaded into your environnment, it should be if you follow the pipeline.
  overall = apply(tmp,1, mean ,na.rm=T)
  
  
  vals_per_genes <- do.call(cbind,vals)
  rownames(vals_per_genes)<-names(feat_indexed_probes_regions)
  foo <- intersect(names(overall),rownames(vals_per_genes))
  vals_per_genes <- cbind(vals_per_genes,overall[foo])
  colnames(vals_per_genes) <-c(binlist,"overall")
  rownames(vals_per_genes) <- foo
  vals_per_genes = vals_per_genes[order(vals_per_genes[,7]),]
  
  
  return(vals_per_genes)
  
}

#################heatmap######################

meth_heatmap <- function(data = means, dendrogram = "none", Rowv = NULL, Colv = NULL, main= "", ...){
  
  
  data = data[,-7]
  
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(20)
  
  foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, trace="none", col=cols, mar=c(10,5), useRaster=TRUE, main= main)
  
  
}



###################preprocess cnv ##############


process_cnv <- function(data_cnv = cnv_lusc$data,
                        treshold = 0.3){
  
  processed_cnv<- epimedtools::monitored_apply(data_cnv,1,
                                               function(gene){
                                                 sapply(gene,
                                                        function(patient){
                                                          if(patient <= abs(treshold) & patient >= -abs(treshold)){patient <- 1}
                                                          else{patient <- 0}
                                                          
                                                        })
                                               })
  return(t(processed_cnv))
}



################## preprocess meth data#################



process_meth_data<-function(cnv_processed_data = cnv_processed,
                            meth_data = meth_lusc$data,
                            probes_index = feat_indexed_probes,
                            features_list = features){
  
  
  
  targeted_cnv <- cnv_processed_data[intersect(rownames(cnv_processed_data),rownames(features_list)),
                                     intersect(colnames(cnv_processed_data),colnames(meth_data))]
  
  methvals_per_genes <- epimedtools::monitored_apply(mod = 10, t(t(names(probes_index))), 1, function(f) {
    
    
    tmp_probes <- intersect(probes_index[[f]],rownames(meth_data))
    selected_meth <- meth_data[tmp_probes,intersect(colnames(targeted_cnv),colnames(meth_data))]  #get data for the f gene
    
    
    
    dummy_matrix <- matrix(nrow = nrow(selected_meth),ncol = ncol(selected_meth),
                           dimnames = list(rownames(selected_meth),colnames(selected_meth)))
    
    
    tmp<-apply(t(t(colnames(selected_meth))),1,function(i) {
      if(targeted_cnv[f,i] != 1){dummy_matrix[,i] <- NA }      #process values per patient for the f gene
      else{dummy_matrix[,i]<-selected_meth[,i]}
    })
    
    
    methvals <- do.call(cbind,tmp)                      #processed meth values for a given f gene
    return(methvals)
    
  })
  
  meth_processed <- do.call(rbind,methvals_per_genes) 
  colnames(meth_processed)<- colnames(targeted_cnv)
  
  return(meth_processed)
}

