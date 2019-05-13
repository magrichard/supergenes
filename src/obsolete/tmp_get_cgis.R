get_cgis_probes <- function(cgi.coordinates = cgi_coordinates,
                               cgi_platform = cgi_pf,
                               meth_platform = meth_lusc$platform){
  

  
  
  # get features type
  
  types<-unique(meth_platform[,"Feature_Type"])

  
  # index data per chromosome to reduce compute time
  
  chrs <- unique(cgi.coordinates[, 1])
  
  print("Indexing CGI per chromosome...")
  chrs_indexed_probes <- lapply(chrs, function(chr) {
    
    idx <- rownames(cgi.coordinates)[cgi.coordinates[["chrs"]] %in% chr]
    ret <- cgi.coordinates[idx, ]
    return(ret)
  })
  names(chrs_indexed_probes) <- chrs
  
  
  
  tmp <- meth_platform[,c("Chromosome","Feature_Type")]
  print("Indexing probes per chromosome")
  chrs_indexed_type <- lapply(chrs, function(chr) {
    
    idx <- rownames(tmp)[tmp[["Chromosome"]] %in% chr]
    ret <- tmp[idx, ]
    return(ret)
  })
  names(chrs_indexed_type) <- chrs
  
  
  cgi_pf <- cgi_platform[which(cgi_platform[,"chr"] != "*"),]
  
  probes_per_cgi<-epimedtools::monitored_apply(t(t(rownames(cgi_pf))),1,mod=50, function(cgi){
    
    # get probes on ith CGI
    
    chr<- cgi_pf[cgi,"chr"]
    probes_on_chr <- chrs_indexed_probes[[as.character(chr)]]
    tmp_probes<-rownames(probes_on_chr[which(probes_on_chr[["cgi_start"]] == cgi_pf[cgi,"start"]),])
    probes_type <- chrs_indexed_type[[as.character(chr)]][tmp_probes,]
    
    # Index probes per feature type (island, shore...)
    
    ret<- apply(t(t(types)),1,function(type){
      
      foo<-assign(paste0(type),rownames(probes_type[which(probes_type[,"Feature_Type"] == type),]))
      
      type_probes <- foo[which(!is.na(foo))]
      
      return(type_probes)
      
    })
    
    names(ret)<- types
    
    return(ret)
    
    
  })
  
  names(probes_per_cgi)<- rownames(cgi_pf)
  
  return(probes_per_cgi)
}
  

    
