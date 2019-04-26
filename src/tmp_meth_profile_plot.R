
plot_gene_cgi <- function(selected_gene = "OTX1",
                           window = c(100000,100000),
                           meth_platform = meth_lusc$platform,
                           meth_data = meth_lusc$data,
                           cgi_platform = cgi_pf,
                           probes_chip_index = chip_index,
                           features_list = features){

  

  TSS <- features_list[selected_gene,"TSS"]
  
  
  ## get closest cgi + cgi in range
  
  cgi_chr <- cgi_platform[which(as.character(cgi_platform[,"chr"]) == features_list[selected_gene,1]),]
  tmp <- cgi_chr[which(abs(cgi_chr[,"center"]-TSS) < mean(window)+5000),]
  cgi_of_interest<- cbind(tmp,abs(tmp[,"center"]-TSS))
  colnames(cgi_of_interest)<-c(colnames(tmp),"dist to TSS")
  closest_cgi <-cgi_chr[which(abs(cgi_chr[,"center"]-TSS)==min(abs(cgi_chr[,"center"]-TSS),na.rm=T)),]
  
  ## get probes
  
  tmp_probes <- dmprocr::get_probe_names(features_list[selected_gene,], meth_platform, "Chromosome", "Start", up_str=100000, dwn_str=100000)
  sub_epic <- meth_platform[tmp_probes, 2:3]
  
  ## cols of probes dots
  
  chip_index_of_interest <- probes_chip_index[intersect(rownames(sub_epic),rownames(probes_chip_index)),]
  
  cols<-sapply(rownames(sub_epic), function(probe){
    
    if(chip_index_of_interest[probe,"epic850k"] == 1)
    {col ="gray55"}
    
    if(chip_index_of_interest[probe,"epic450k"] == 1)
    {col = "darkolivegreen"}
    
    
    if(chip_index_of_interest[probe,"epic27k"] == 1)
    {col = "darkorchid"}
    
    
    if(chip_index_of_interest[probe,"epic850k"] == 1 & chip_index_of_interest[probe,"epic450k"] == 1)
    {col = "darkcyan"}
    
    
    if(chip_index_of_interest[probe,"epic27k"] == 1 & chip_index_of_interest[probe,"epic850k"] == 1)
    {col = "darkred"}
    
    
    if(chip_index_of_interest[probe,"epic27k"] == 1 & chip_index_of_interest[probe,"epic450k"] == 1)
    {col = "goldenrod"}
    
    
    if(chip_index_of_interest[probe,"epic27k"] == 1 & chip_index_of_interest[probe,"epic450k"] == 1 & chip_index_of_interest[probe,"epic850k"] == 1)
    {col = "black"}
    
    
    return(col)
  })
  
  layout(matrix(1:2,1), respect=TRUE)

  plot(sub_epic$Start, rep(1,length(sub_epic$Start)), pch=19, xlim=c(TSS-window[1],TSS+window[2]), cex=1.1, yaxt="n",
       main = paste0(selected_gene,": Probes repartition among CGI over selected window \n","nprobes = ",nrow(sub_epic),", ncgi = ",nrow(cgi_of_interest)),
       xlab="Coordinates (in bp)",
       ylab="",
       ylim=c(0,2),
       col = cols)
  
  
  text(TSS,0.4,labels = paste0("TSS : ", TSS),cex = 0.7)
  points(TSS, 0.5, pch=9, col="red")
  abline(h=0.5,lty=3)

  if(nrow(cgi_of_interest)>0){
    apply(t(t(rownames(cgi_of_interest))),1,function(cgi){
      
      rect(cgi_of_interest[cgi,"start"],1.4,cgi_of_interest[cgi,"end"],1.2,
           col = "red",
           border="black")
      
    })
  }
  
  ## Indexing probes per cgi
  
  tmp_cgi_index <- apply(t(t(rownames(sub_epic))),1, function(probe){
                foo <-apply(t(t(rownames(cgi_of_interest))),1,function(cgi){
                    
                        if(sub_epic[probe,"Start"] < cgi_of_interest[cgi,"end"] & sub_epic[probe,"Start"] >= cgi_of_interest[cgi,"start"]){
                         rownames(cgi_of_interest[cgi,])}
                  })
                
                                  index<-unlist(foo)
                                  if(!is.null(index)){
                                  names(index)<- probe
                                  }
                                  return(index)
  })
  
  cgi_index <- sort(unlist(tmp_cgi_index))
  
  methvals_of_interest <- meth_data[intersect(rownames(meth_data),names(cgi_index)),]
  methvals_of_interest<- methvals_of_interest[names(cgi_index),]
  
  
  colors=c("cyan", "black", "red")
  cols = colorRampPalette(colors)(100)
  
  breaks <- c(seq(0,0.33,length=35),seq(0.34,0.66,length=33),seq(0.67,1,length=33))
  image(methvals_of_interest, axes=FALSE, col=cols, main="Methylation",breaks = breaks)
  
}
