
  

get_indexed_binreg<-function(features_list = features,
                             epimed = epic,
                             platform_regions = platform){
  
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


  

  print(noquote(paste0("Indexing probes per bins...")))



  tmp_bins <- epimedtools::monitored_apply(features_list, 1, function(gene) {


    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[as.character(chr)]]
    tmp_probes <- feat_indexed_probes[[gene[[4]]]]
    sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname,"end")]
    
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


  
  regions_name <-unique(as.character(platform_regions[,"genomic_feature"]))[unique(as.character(platform_regions[,"genomic_feature"]))!= "P2000"]
  
  regions <- epimedtools::monitored_apply(t(t(rownames(features_list))),mod=20,1, function(gene){
    print(noquote(paste0("Indexing probes per regions for ",gene,"...")))
    
    
    TSS <- features_list[gene,"TSS"]
    
    
    tmp_probes <- feat_indexed_probes[[gene]]
    tmp_pf <- platform_regions[which(platform_regions[,"gene"]==gene),]
    sub_epic <- epimed[intersect(tmp_probes,rownames(data)),c("start","end")]
    
    if(features_list[gene,"strand"] == "+"){
      pf_gene <- tmp_pf[which(tmp_pf[,"end"] > TSS),]
    }else{pf_gene<- tmp_pf[which(tmp_pf[,"end"]< TSS ),]}
    
    
    
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
    tmp <- unlist(c(bins[gene],regions[gene]),recursive=FALSE,use.names=FALSE)
    names(tmp)<-c(names(bins[[gene]]),names(regions[[gene]]))
    return(tmp)
  })
  names(feat_indexed_probes_binreg)<- rownames(features_list)
  
  return(feat_indexed_probes_binreg)
  
}
