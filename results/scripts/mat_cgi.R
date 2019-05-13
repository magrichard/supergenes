
get_cgi_mat<-function(features_list = features,
                      window = c(5000,5000),
                      bwidth = 25,
                      cgi_platform = cgi_pf){

  nbins <- sum(window)/bwidth
  bins_coordinates <- seq(-abs(window[1]),window[2],bwidth)

  
  cgi_mat <- matrix(NA,nrow = nrow(features), ncol = nbins)
  rownames(cgi_mat)<- rownames(features_list)
  
  
  chrs <- unique(cgi_platform[, 1])
  
  print("Indexing CGI per chromosome...")
  
  chrs_indexed_cgi <- lapply(chrs, function(chr) {
    
    idx <- rownames(cgi_platform)[cgi_platform[["chr"]] %in% chr]
    ret <- cgi_platform[idx, ]
    return(ret)
  })
  names(chrs_indexed_cgi) <- chrs

  #gene = "OTX1"
  
  for(gene in 1:nrow(features)){
  
  
    bins_coordinates <- seq(-abs(window[1]),window[2],bwidth)
  
    if(features[gene,"strand"]=="-"){bins_coordinates <- sort(bins_coordinates,decreasing = TRUE)}
    
  TSS <- features_list[gene,"TSS"]
  
  
  cgi_on_chr <- chrs_indexed_cgi[[as.character(features[[gene,1]])]]
  
  tmp_cgi_of_interest <- cgi_on_chr[which(abs(cgi_on_chr[,"center"]-TSS) < mean(window)+3500),]
  cgi_of_interest<- cbind(tmp_cgi_of_interest,abs(tmp_cgi_of_interest[,"center"]-TSS))
  colnames(cgi_of_interest)<-c(colnames(tmp_cgi_of_interest),"dist to TSS")
  if(nrow(cgi_of_interest)< 1){next}
  
    for(cgi in 1:nrow(cgi_of_interest)){
      
     for (i in 2:length(bins_coordinates)){
        
        if((cgi_of_interest[cgi,"start"] <= TSS + bins_coordinates[i-1]) & (TSS + bins_coordinates[i-1] <= cgi_of_interest[cgi,"end"])
           |(cgi_of_interest[cgi,"start"] >= TSS + bins_coordinates[i-1]) & (TSS + bins_coordinates[i-1] >= cgi_of_interest[cgi,"end"]))
              {cgi_mat[gene,i-1]<- 1}
        }
    }
  }
          
  cgi_mat[which(is.na(cgi_mat))]<- 0
  
  d = dist(cgi_mat)
  hc_row = hclust(d, method="complete")
  Rowv = as.dendrogram(hc_row)
  dendrogram="row"
  
  
  meth_heatmap(cgi_mat,main = paste0("CGI Binary Matrix : ",window[1],";",window[2]),
               dendrogram = dendrogram,
               Rowv=Rowv)
  
  return(cgi_mat)
  
}
  
  
  



  
  
  
  
  
  