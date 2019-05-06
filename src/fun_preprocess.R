###################################  Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################


get_features <- function(targeted_genes, study, up_str = 2500, dwn_str = 2500, nb_probe_min = 1) { 
  print(noquote("Creating a list of features..."))
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

  print(noquote("Indexing probe by features"))
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
  gc()

  return(sub_features)
}




############################### get_features_bins : same as prior + indexing probes per bins #####################

get_features_bins <- function(targeted_genes, study = trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1, ...) { #### In my uses-case, 1 for LUAD & 2 for LUSC
  print("Creating a list of features...")
  features <- study$platform[intersect(rownames(study$platform), targeted_genes), ]
  
  
  
  
  print(noquote("Indexing probes by features and by bins..."))
  
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
  
  gc()
  
  return(feat_indexed_probes_bin)
}





################# Index regions per features################
get_features_regions <- function(features_list = features,
                                 platform_regions = platform,
                                 data = meth_lusc$data,
                                 epimed_database = epic,
                                 indexed_probes = feat_indexed_probes ){
  
  regions_name <-unique(as.character(platform[,"genomic_feature"]))
  
  feat_indexed_probes_regions <- lapply(1:nrow(features_list), function(index){
    
    
    gene <- names(feat_indexed_probes[index])
    tmp_probes <- feat_indexed_probes[[gene]]
   
    sub_epic <- epic[intersect(tmp_probes,rownames(data)),c("start","end")]
    
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

get_feat_indexed_DMRs <- function(DMRtable = DMR_table,
                                  features_list = features,
                                  regions_coordinates = platform,
                                  treshold = 100000){
  
  DMRtable<- DMRtable[order(DMRtable[,"DMR_id"]),]
  
  print("Extractig regions coordinates per available genes...")
  
  pf<- lapply(rownames(features_list), function(gene){
    regions_coordinates[which(regions_coordinates[,"gene"]==gene),]
  })
  names(pf)<- rownames(features_list)
  
  
  
  print("Indexing DMRs per genes...")
  
  feat_indexed_DMRs <- epimedtools::monitored_apply(t(t(rownames(features_list))),mod=10,1, function(gene){
    
    print(noquote(paste0("Indexing DMRs for ", gene ,"...")))
    
    DMRs_of_interest <- DMRtable[which(DMRtable[,1]==features[gene,1]),]
    
    DMR_in_range <- lapply(1:nrow(DMRs_of_interest),function(index){
      DMR <- DMRs_of_interest[index,]
      
      if(nrow(DMRs_of_interest) > 0){
        if(as.numeric(DMR[["start"]]) > features_list[gene,"TSS"]){
          foo <- abs(as.numeric(DMR[["start"]])-features_list[gene,"TSS"])
          return(foo)
        }
        
        if(as.numeric(DMR[["start"]]) < features_list[gene,"TSS"]){
          foo <- abs(features_list[gene,"TSS"]-as.numeric(DMR[["end"]]))
          return(foo)
        }
      }
      
    })
    if (nrow(DMRs_of_interest) > 0){
      names(DMR_in_range)<-DMRs_of_interest[["DMR_id"]]
      DMR_candidates <- unlist(DMR_in_range)[DMR_in_range <= treshold]
    }
    else{DMR_candidates<- NULL}
    
    return(DMR_candidates)
    
  })
  
  names(feat_indexed_DMRs)<- rownames(features_list)
  
  
  return(feat_indexed_DMRs)
  
  gc()
}

######################## get_meth_diff ##############
get_differential_values <- function(meth_tumoral,meth_healthy){
  
  mean_values_normal_tissues <- apply(meth_healthy,1,mean,na.rm=T)
  meth_diff <- t(sapply(1:nrow(meth_tumoral), function(i) meth_tumoral[i,] - mean_values_normal_tissues[i]))
  rownames(meth_diff)<-rownames(meth_tumoral)
  return(meth_diff)
}



############catalog##################

catalog <- function(selected_features = features) {
  available_genes <- list(catalog = "catalog", n_genes = "n_genes")
  
  available_genes[["catalog"]] <- noquote(rownames(selected_features))
  available_genes [["n_genes"]] <- noquote(paste("Genes selected :", length(available_genes[[1]]), sep = " "))
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


reduce_map <- function(map_to_reduce = feat_indexed_probes_bin, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"))  {
  
  bin_indexed_probes = lapply(map_to_reduce, function(g){
    g[binlist]
  })
  bin_indexed_probes = unlist(bin_indexed_probes, recursive = FALSE)
  
}

####################### get means per bins per genes #############################

#means_per_bins_per_genes_per_patient = (meth_lusc$data,binmap,mean,na.rm=T)

subset_vals_per_bins<-function(data = meth_lusc$data,
                               values_per_patient = means_per_regions_per_genes_per_patient,
                               binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"),
                               fun = mean,
                               probes_index = feat_indexed_probes, ...){
  
  names<-names(which(unlist(lapply(probes_index,is.null))==FALSE))
  
  
  values_per_bins <- apply(values_per_patient,1,mean,na.rm=T,...) 
  ret <- list()
  
  vals <- sapply(binlist, function(bin){
    ret[[bin]] <- values_per_bins[endsWith(names(values_per_bins),bin)]    
    return(ret)
    
  }
  )
  
  
  tmp = reduce_rows(data, feat_indexed_probes, fun, na.rm=T,...) ### note that you need feat_indexed_probes loaded into your environnment, it should be if you follow the pipeline.
  overall = apply(tmp,1, mean , na.rm=T,...)
  
  
  vals_per_genes <- do.call(cbind,vals)
  rownames(vals_per_genes)<-names
  foo <- intersect(names(overall),rownames(vals_per_genes))
  vals_per_genes <- cbind(vals_per_genes,overall[foo])
  colnames(vals_per_genes) <-c(binlist,"overall")
  if(!is.null(foo)){
  rownames(vals_per_genes) <- foo}
  
  
  return(vals_per_genes)
  
  gc()
  
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
  
  gc()
}



################## preprocess meth data#################



process_meth_data<-function(cnv_processed_data = cnv_processed,
                            meth_data = meth_lusc$data,
                            probes_index = feat_indexed_probes,
                            features_list = features){
  
  
  
  targeted_cnv <- cnv_processed_data[intersect(rownames(cnv_processed_data),rownames(features_list)),
                                     intersect(colnames(cnv_processed_data),colnames(meth_data))]
  
  methvals_per_genes <- epimedtools::monitored_apply(mod = 10, t(t(rownames(features_list))), 1, function(f) {
    
    
    tmp_probes <- intersect(probes_index[[f]],rownames(meth_data))
    
    if(length(tmp_probes >= 1)){ 
      
      selected_meth <- meth_data[tmp_probes,intersect(colnames(targeted_cnv),colnames(meth_data))]  #get data for the f gene
      
      
      if (class(selected_meth) == "numeric"){
        selected_meth <- t(selected_meth)
      }
      
      
      
      dummy_matrix <- matrix(nrow = nrow(selected_meth),ncol = ncol(selected_meth),
                             dimnames = list(rownames(selected_meth),colnames(selected_meth)))
      
      
      tmp <- apply(t(t(colnames(selected_meth))),1,function(i) {
        if(targeted_cnv[f,i] != 1){
          dummy_matrix[,i] <- NA 
        }      #process values per patient for the f gene
        else{
          dummy_matrix[,i]<-selected_meth[,i]
        }
      })
      
      if(class(tmp) == "list"){
        methvals <- do.call(cbind,tmp)                      #processed meth values for a given f gene
      }
      
      else{methvals <- tmp}
      
      return(methvals)
      
    }
    
  })
  
  meth_processed <- do.call(rbind,methvals_per_genes) 
  colnames(meth_processed)<- colnames(targeted_cnv)
  meth_processed[unique(rownames(meth_processed)),]
  return(meth_processed)
  
  gc()
}



##################### get probes indexed by bins and by regions ############



get_indexed_binreg<-function(features_list = features,
                             epimed = epic,
                             platform_regions = platform,
                             data = meth_lusc$data){
  
  pf_chr_colname <- "seqnames"
  pf_pos_colname <- "start"
  chrs <- unique(features_list[, 1])
  chrs_indexed_epic <- lapply(chrs, function(chr) {
    idx <- rownames(epimed)[epimed[[pf_chr_colname]] %in% chr]
    ret <- epimed[idx, ]
    return(ret)
  })
  names(chrs_indexed_epic) <- chrs
  
  
  
  print(noquote(paste0("Indexing probes per bins...")))
  
  
  
  tmp_bins <- epimedtools::monitored_apply(features_list, 1, function(gene) {
    
    
    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[chr]]
    tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, 7500, 7500)
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
    
    
    ret <- list(bin1 = bin1,
                bin2 = bin2,
                bin3 = bin3)
    
    return(ret)
    
  })
  
  bins <- apply(t(t(names(tmp_bins))),1,function(gene){
    lapply(tmp_bins[[gene]],unlist,recursive = FALSE)
  })
  names(bins)<-rownames(features_list)
  
  
  
  regions_name <- c("INTER","CDS","INTRON", "UTR5","UTR3")
  
  regions <- epimedtools::monitored_apply(t(t(rownames(features_list))),mod=20,1, function(gene){
    print(noquote(paste0("Indexing probes per regions for ",gene,"...")))
    
    
    TSS <- features_list[gene,"TSS"]
    
    
    tmp_probes <- feat_indexed_probes[[gene]]
    tmp_pf <- platform_regions[which(platform_regions[,"gene"]==gene),]
    sub_epic <- epimed[intersect(tmp_probes,rownames(data)),c("start","end")]
    
    if(features_list[gene,"strand"] == "+"){
      pf_gene <- tmp_pf[which(tmp_pf[,"end"] > TSS),]
    }
    
    if(features_list[gene,"strand"] == "-"){
      pf_gene <- tmp_pf[which(tmp_pf[,"start"] < TSS),]
    }
    
    
    
    pf_regions = lapply(regions_name, function(region){
      regions<-pf_gene[which(pf_gene[,"genomic_feature"]==region),]
      return(regions)
    })
    
    names(pf_regions) <- regions_name
    
    
    
    
    INTER = vector()
    P2000 = vector()
    UTR5 = vector()
    INTRON = vector()
    CDS = vector()
    UTR3 = vector()
    
    for (i in 1:nrow(sub_epic)){
      
      if (nrow(pf_regions$INTER)>0)
        if(sum(apply(pf_regions$INTER,1, function(window){
          sub_epic[i,"start"] > as.numeric(window[["start"]]) & sub_epic[i,"end"] < as.numeric(window[["end"]]) }),na.rm=T)>=1)
        {INTER[i]<-rownames(sub_epic[i,])}
      
      
      
      if (nrow(pf_regions$UTR5)>0)
        if(sum(apply(pf_regions$UTR5,1, function(window){
          sub_epic[i,"start"] > as.numeric(window[["start"]]) & sub_epic[i,"end"] < as.numeric(window[["end"]]) }),na.rm=T)>=1)
        {UTR5[i]<-rownames(sub_epic[i,])}
      
      
      
      if (nrow(pf_regions$INTRON)>0)
        if(sum(apply(pf_regions$INTRON,1, function(window){
          sub_epic[i,"start"] > as.numeric(window[["start"]]) & sub_epic[i,"end"] < as.numeric(window[["end"]]) }),na.rm=T)>=1)
        {INTRON[i]<-rownames(sub_epic[i,])}
      
      
      if (nrow(pf_regions$CDS)>0)
        if(sum(apply(pf_regions$CDS,1, function(window){
          sub_epic[i,"start"] > as.numeric(window[["start"]]) & sub_epic[i,"end"] < as.numeric(window[["end"]]) }),na.rm=T)>=1) 
        {CDS[i]<-rownames(sub_epic[i,])}
      
      if (nrow(pf_regions$UTR3)>0)
        if(sum(apply(pf_regions$UTR3,1, function(window){
          sub_epic[i,"start"] > as.numeric(window[["start"]]) & sub_epic[i,"end"] < as.numeric(window[["end"]]) }),na.rm=T)>=1)
        {UTR3[i]<-rownames(sub_epic[i,])}
      
      INTER <- INTER[!is.na(INTER)]
      INTRON <- INTRON[!is.na(INTRON)]
      CDS <- CDS[!is.na(CDS)]
      UTR3 <- UTR3[!is.na(UTR3)]
      UTR5 <- UTR5[!is.na(UTR5)]
      
      
    }
    
    ret <- list(sub_epic = sub_epic, INTER = INTER, UTR5 = UTR5, CDS = CDS, INTRON = INTRON, UTR3 = UTR3)
    return(ret)
    
  })
  names(regions)<- rownames(features_list)
  
  
  feat_indexed_probes_binreg <- apply(t(t(rownames(features_list))),1,function(gene){
    tmp <- unlist(c(regions[gene],bins[gene]),recursive=FALSE,use.names=FALSE)
    names(tmp)<-c(names(regions[[gene]]),names(bins[[gene]]))
    return(tmp)
  })
  names(feat_indexed_probes_binreg)<- rownames(features_list)
  
  return(feat_indexed_probes_binreg)
  
}


################# Get probes indexed by cgi features ###################

get_indexed_cgi <- function(window = c(100000,100000),
                            features_list = features,
                            cgi_platform = cgi_pf,
                            meth_data = meth_lusc$data,
                            meth_platform = meth_lusc$platform,
                            genes_platform = trscr_lusc$platform){
  
  
  cgi_indexed_probes<- epimedtools::monitored_apply(t(t(rownames(features_list))),mod = 20,1,function(gene){
    print(noquote(paste0("Indexing probes for gene ",gene," ...")))
    
    chr<-genes_platform[gene,1]
    TSS <- features_list[gene,"TSS"]
    
    meth_platform_chr<- meth_platform[rownames(meth_platform)[meth_platform[[1]] %in% chr],]
    tmp_probes <- dmprocr::get_probe_names(features_list[gene,], meth_platform_chr, "Chromosome", "Start", window[1], window[2])
    
    ## get probe for gene
    
    sub_platform <- meth_platform_chr[tmp_probes, c("Chromosome","Start","Feature_Type")]
    
    
    
    for (i in which(sub_platform[,3]==".")){
      if(sub_platform[i,"Start"] > TSS){sub_platform[i,"Feature_Type"] <- "upper_opensea"}
      if(sub_platform[i,"Start"] < TSS){sub_platform[i,"Feature_Type"] <- "lower_opensea"}  
      ### encode "." as opensea
    }
    features_type_names <- unique(sub_platform[,"Feature_Type"])
    
    foo<-apply(t(t(features_type_names)),1,function(type){
      
      
      tmp_type<-assign(paste0(type),rownames(sub_platform[which(sub_platform[,"Feature_Type"]==type),]))
      if(!is.list(tmp_type)){assign(paste0(type),list(tmp_type))}
      else{assign(paste0(type),tmp_type)}
      ## create a list of probes for each level of "Feature_Type" + very very very smart trick to deal with 1 length list becoming
      ## character
      

      
    })
    
    for (i in 1:length(foo)){
      if(is.list(foo[[i]])){foo[i]<-foo[[i]]}
    }
      
      
      
    
    names(foo)<- features_type_names
    
    
    
    
    
    ret = list(sub_epic = sub_platform,
               lower_opensea = foo$lower_opensea,
               S_shore = foo$S_Shore,
               S_shelf = foo$S_Shelf,
               Island = foo$Island,
               N_shelf = foo$N_Shelf,
               N_shore = foo$N_Shore,
               upper_opensea = foo$upper_opensea)
    
    return(ret)
    
    
  })
  
  names(cgi_indexed_probes)<- rownames(features_list)
  return(cgi_indexed_probes)
  
}

############# Index probes per DMR #################

get_dmrs_indexed_probes <-function(DMRtable = DMR_table,
                                   meth_platform = meth_lusc$platform){
  
  
  
  
  chrs <- unique(DMRtable[, 1])
  chrs_indexed_probes<- lapply(chrs, function(chr) {
    print(paste0("Indexing probes per chromosome :", as.character(chr),"..."))
    idx <- rownames(meth_platform)[meth_platform[[1]] %in% chr]
    ret <- meth_platform[idx, ]
    return(ret)
  })
  names(chrs_indexed_probes) <- chrs
  
  
  print("Indexing probes per DMRs...")
  probes_per_dmr<-epimedtools::monitored_apply(t(t(rownames(DMRtable))),mod=500,1, function(i){
    
    DMR<-DMRtable[i,]
    probes_on_chr<- chrs_indexed_probes[[DMR[["chr"]]]]
    
    
    selected_probes <- rownames(probes_on_chr[which(probes_on_chr[,2] >= DMR[["start"]] & probes_on_chr[,3] <= DMR[["end"]]),
                                              ])
    return(selected_probes)
    
  })
  
  names(probes_per_dmr)<- rownames(DMRtable)
  
  return(probes_per_dmr)
  
}
