
get_features_regions <- function(platform_regions = platform, epimed_database = epic, indexed_probes = feat_indexed_probes ){

regions_name <-unique(as.character(platform[,"genomic_feature"]))

  feat_indexed_probes_regions <- epimedtools::monitored_apply(features[1:50,], 1, function(gene) {


    tmp_probes <- unlist(indexed_probes[gene])
  
    sub_epic <- epic[tmp_probes,c("start","end")]
    pf_gene <- platform[which(platform[,"gene"]==gene),]


    pf_regions = lapply(regions_name, function(region){
  
      regions<-pf_gene[which(pf_gene[,"genomic_feature"]==region),]
      return(regions)
    })

    names(pf_regions) <- regions_name

  per_regions_type_coordinates <-sapply(pf_regions, function(pf){
    if (nrow(pf)>1){
      per_regions_coordinates = apply(pf,1,function(region){
        range <- range(region["start"],region["end"])
        seq(range[1],range[2])
        
      })
    }
  })


  INTER = list()
  P2000 = list()
  UTR5 = list()
  INTRON = list()
  CDS = list()
  UTR3 = list()

  for (i in 1:nrow(sub_epic)){

    if (!is.null(per_regions_type_coordinates$INTER)){
      if(sum(sapply(per_regions_type_coordinates$INTER, function(f) is.element(sub_epic[i,"start"],f)))>=1)   
      {INTER[[i]]<-rownames(sub_epic[i,])}
    }
    
    if (!is.null(per_regions_type_coordinates$P2000))
      if(sum(sapply(per_regions_type_coordinates$P2000, function(f) is.element(sub_epic[i,"start"],f)))>=1)   
      {P2000[[i]]<-rownames(sub_epic[i,])}
    
    if (!is.null(per_regions_type_coordinates$CDS))
      if(sum(sapply(per_regions_type_coordinates$CDS, function(f) is.element(sub_epic[i,"start"],f)))>=1)   
      {CDS[[i]]<-rownames(sub_epic[i,])}
    
    if (!is.null(per_regions_type_coordinates$INTRON))
      if(sum(sapply(per_regions_type_coordinates$INTRON, function(f) is.element(sub_epic[i,"start"],f)))>=1)   
      {INTRON[[i]]<-rownames(sub_epic[i,])}
    
    if (!is.null(per_regions_type_coordinates$`5UTR`))
      if(sum(sapply(per_regions_type_coordinates$`5UTR`, function(f) is.element(sub_epic[i,"start"],f)))>=1)   
      {UTR5[[i]]<-rownames(sub_epic[i,])}
    
    if (!is.null(per_regions_type_coordinates$`3UTR`))
      if(sum(sapply(per_regions_type_coordinates$`3UTR`, function(f) is.element(sub_epic[i,"start"],f)))>=1)   
      {UTR3[[i]]<-rownames(sub_epic[i,])}
  
  }

  ret <- list(sub_epic = sub_epic, whole_feature = tmp_probes, INTER = INTER, P2000 = P2000, UTR5 = UTR5, CDS = CDS, INTRON = INTRON, UTR3 = UTR3)
  return(ret)
  })


}
  
  
  
  