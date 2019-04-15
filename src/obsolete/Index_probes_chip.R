probes_info<-matrix(nrow=length(unique(c(rownames(epic),rownames(epic27k),rownames(epic450k)))), ncol=3)
rownames(probes_info) <- unique(c(rownames(epic),rownames(epic27k),rownames(epic450k)))
colnames(probes_info)<- c("epic850k","epic450k","epic27k")

for (i in 1:nrow(probes_info)){
  
  probe<-rownames(probes_info)[i]
  
  if (probe %in% rownames(epic)){
    probes_info[i,1]<- 1
  }
  else {probes_info[i,1]<- 0}
  
  if (probe %in% rownames(epic450k)){
    probes_info[i,2]<- 1
  }
  else {probes_info[i,2]<- 0}
  
  if (probe %in% rownames(epic27k)){
    probes_info[i,3]<- 1
  }
  else {probes_info[i,3]<- 0}

}

probes_info <- data.frame(probes_info)
write.csv(probes_info, file ="probes_chip_index")
