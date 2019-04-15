cgis<-meth_lusc$platform[,"CGI_Coordinate"]
coordinates <- sapply(1:nrow(meth_lusc$platform),function(x){
  str_split_fixed(gsub("^([^:]*:[^:]*):","",cgis[x]), "-", 2)
  })
cgi_coordinates<-t(coordinates)


cgi <- cbind(as.numeric(cgi_coordinates[,1]),as.numeric(cgi_coordinates[,2]))
rownames(cgi)<-rownames(meth_lusc$data)
colnames(cgi)<- c("cgi_start","cgi_end")
head(cgi)
