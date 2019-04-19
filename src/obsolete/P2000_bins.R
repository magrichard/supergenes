
tmp <- lapply(rownames(features), function(gene){
  platform[which(platform[,"gene"]==gene),]              ### Get informations for targeted_genes
})
pf_tmp <- do.call(rbind,tmp)

P2000 <- pf_tmp[which(pf_tmp[,"genomic_feature"]=="P2000"),]
rownames(P2000) <- P2000[,"gene"] 
pf <- pf_tmp[which(pf_tmp[,"genomic_feature"]!="P2000"),] 

nbins = 5
binswidth = 2000/nbins

binsnames <- c()
for (i in 1:nbins){
  binsnames[i]<- paste0("P2000.bin",i)   # create bin names
}




  
  
  
  end<-P2000[gene,"end"]
  start <- P2000[gene,"start"]
  
  
  bins_ends <- sort(seq(end,start+(binswidth-1),-binswidth))
  bins_starts <- seq(start,end-(binswidth-1),binswidth)
  
  chr <- rep(P2000[gene,"chromosome"],nbins)
  gene_name <- rep(P2000[gene,"gene"],nbins) 
  
  foo <- data.frame(chromosome = chr,
                    start = bins_starts,
                    end = bins_ends,
                    genomic_feature = binsnames,
                    gene = gene_name)
  
    return(foo)

z<- c()

  for (i in 1:nrow(P2000)){
    z[i]<- P2000[i,"start"]-P2000[i,"end"]
    
  }








