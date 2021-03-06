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

#DMRs_tmp2 <- get_indexed_DMRs(DMRtable = DMR_table)

