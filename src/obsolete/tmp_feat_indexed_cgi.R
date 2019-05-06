source("~/projects/supergenes/src/import_data.R")


get_indexed_cgi <- function(window = c(100000,100000),features_list = features,cgi_platform = cgi_pf,meth_data = meth_lusc$data,meth_platform = meth_lusc$platform,gene_platform = trscr_lusc$platform){
  
  
  
  cgi_indexed_probes<- epimedtools::monitored_apply(t(t(rownames(features_list))),mod = 20,1,function(gene){
    print(noquote(paste0("Indexing probes for gene ",gene," ...")))
    
    chr<-gene_platform[gene,1]
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
      assign(paste0(type),rownames(sub_platform[which(sub_platform[,"Feature_Type"]==type),])) 
      ## create a list of probes for each level of "Feature_Type"
      
    })
    
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













