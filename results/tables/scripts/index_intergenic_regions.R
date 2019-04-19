pf_inter_start <- platform[which(platform[,"genomic_feature"]=="INTER"),]
pf_inter_end <- platform[which(platform[,"genomic_feature"]=="INTER"),]

pf_inter_duplicated <- epimedtools::monitored_apply(t(pf_inter),1,function(x) rbind(x,x))
pf_inter_duplicated <- data.frame(chromosome = as.factor(pf_inter_duplicated[,1]),
           start = as.numeric(pf_inter_duplicated[,2]),
           end = as.numeric(pf_inter_duplicated[,3]),
           genomic_feature = as.factor(pf_inter_duplicated[,4]),
           gene = as.factor(pf_inter_duplicated[,5]))


inter_lower <- lapply(1:nrow(pf_inter_end),function(index){
  
  inter <- data.frame(pf_inter_end[index,1:4],platform[as.character(as.numeric(rownames(pf_inter_end[index,]))+1),"gene"])
  
    return(inter)
  }
)

inter_lower_dt <- do.call(rbind,inter_lower)
colnames(inter_lower_dt)<- colnames(platform)
inter_lower_dt <- inter_lower_dt[-nrow(inter_lower_dt),]



inter_upper <- lapply(1:nrow(pf_inter_end),function(index){
  
  inter <- data.frame(pf_inter_end[index,1:4],platform[as.character(as.numeric(rownames(pf_inter_end[index,]))-1),"gene"])
  
  return(inter)
  }
)

inter_upper_dt <- do.call(rbind,inter_upper)
colnames(inter_upper_dt)<- colnames(platform)
inter_upper_dt <- inter_upper_dt[-1,]

platform_ <- rbind(inter_upper_dt,inter_lower_dt,platform[which(platform[,"genomic_feature"]!="INTER"),])

apply(platform_,2,class)

write.csv(platform_,"biological_regions.csv")

