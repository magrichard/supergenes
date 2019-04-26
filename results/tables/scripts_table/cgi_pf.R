
cgi_coordinates<- unique(cgi_coordinates)
colnames(cgi_coordinates)<-c("chr","start","end")


mr <- epimedtools::monitored_apply(t(t(rownames(cgi_coordinates))),1,function(cgi){
  
  means<-as.integer(mean(c(cgi_coordinates[cgi,"start"],cgi_coordinates[cgi,"end"])))
  ranges <- cgi_coordinates[cgi,"end"]-cgi_coordinates[cgi,"start"]
  return(c(means,ranges))
})

mr<-data.frame(t(mr))
rownames(mr)<- rownames(cgi_coordinates)


cgi_pf <- cbind(cgi_coordinates,mr)
colnames(cgi_pf)<-c(colnames(cgi_coordinates),"center","width")

write.csv(cgi_pf,"cgi_pf.csv")
