plot_selected_regions<- function(selected_gene = "ANLN",
                                 features_list = features,
                                 DMRtable = DMR_table,
                                 cgi_table = cgi_coordinates,
                                 DMR_map,
                                 pf = platform,
                                 epimed = epic,
                                 chips_index = chip_index,
                                 probes_index = feat_indexed_probes,
                                 window=c(10000,10000)){
  
  ## get nearest DMR
  DMRs_of_interest<-DMRtable[which(DMRtable[,"closest_gene_name"]==selected_gene),]
  #DMRs_of_interest<- DMR_table[intersect(DMR_table[,"DMR_id"],names(DMR_map[[selected_gene]])),]
  
  
  
  ## get data
  
  
  upstream = window[2]
  downstream = window[1]
  pf_gene <- pf[which(pf[,"gene"]==selected_gene),2:4]
  tmp_probes <- feat_indexed_probes[[selected_gene]]
  sub_epic_tmp <- epimed[tmp_probes,c("start","end")]
  sub_epic_tmp <- sub_epic_tmp[intersect(rownames(sub_epic_tmp), rownames(meth_lusc$data)),]
  sub_epic <- cbind(sub_epic_tmp,cgi_table[intersect(rownames(sub_epic_tmp),rownames(cgi_table)),])
  TSS <- features_list[selected_gene,"TSS"]
  chip_index_of_interest <- chips_index[intersect(rownames(sub_epic),rownames(chips_index)),]
  
  
  ## get subset of coordinates
  
  introns_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTRON"),]
  exons_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="CDS"),]
  utr3_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="3UTR"),]
  utr5_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="5UTR"),]
  intergenic_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTER"),]
  p2000_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="P2000"),]
  
  ## get colors
  
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
  
  ## plot probes + TSS
  
  plot(sub_epic$start, rep(3,length(sub_epic$start)), pch=19, xlim=c(TSS-downstream,TSS+upstream), cex=1.2, yaxt="n",
       main = paste0(selected_gene,": Probes repartition among regions (nprobes = ",nrow(sub_epic),")"),
       xlab="Coordinates (in bp)",
       ylab="",
       ylim=c(0,5),
       col = cols)
  
  legend("left", c("P2000", "Exons", "Introns","3'UTR","5'UTR","Inter","hypoDMR","hyperDMR","CGI"),
         inset=c(-0.12,0),
         xpd = TRUE,
         pch=15,
         col=c("grey","darkslategray1","firebrick","darkseagreen","orange","black","blue","red","gray1"),
         bty="n",title = "Region")
  
  legend("topright",
         c("850k","450k","27k","850-450","850-27","450-27","all"),
         inset=c(-0.05,0),
         xpd=TRUE,
         pch = 19,
         col = c("gray55","darkolivegreen","darkorchid","darkcyan","darkred","goldenrod","black"),
         title ="Chip index")
  

  #text(sub_epic$start,rep(4,length(sub_epic$start)), labels = rownames(sub_epic) ,cex = 0.8,srt = 80) ##probes label
  text(TSS,0.2,labels = paste0("TSS : ", TSS),cex = 0.7)
  
  points(TSS, 0.5, pch=9, col="red")
  abline(h=0.5, lty=3, lwd = 1)
  
  
  
  ##plot CGIs
  
 for (i in 1:nrow(sub_epic)){
    
   rect(as.numeric(as.character(sub_epic[i,"cgi_start"])),2.7, as.numeric(as.character(sub_epic[i,"cgi_end"])),2.6,
         col ="gray1")
  }
  
  
  
  
  ## plot regions
  
  
  
  if(nrow(p2000_coordinates) > 0) {
    apply(p2000_coordinates,1, function(p2000){
      rect(p2000[["start"]],1.5, p2000[["end"]],2.5,
           col="grey",
           border ="black")
    })
  }
  
  
  
  if (nrow(exons_coordinates)>0){
    apply(exons_coordinates,1, function(exon){
      rect(exon[["start"]],1.5, exon[["end"]], 2.5,
           col="darkslategray1",
           border="black")
    })
  }
  
  
  if(nrow(introns_coordinates)>0){
    apply(introns_coordinates,1, function(intron){
      rect(intron[["start"]],1.5, intron[["end"]],2.5,
           col="firebrick",
           border="black")
    })
  }
  
  if(nrow(utr3_coordinates)>0){
    apply(utr3_coordinates,1, function(utr3){
      rect(utr3[["start"]],1.5, utr3[["end"]],2.5,
           col="darkseagreen",
           border="black")
      
    })
  }
  
  
  
  if(nrow(utr5_coordinates)>0){
    apply(utr5_coordinates,1, function(utr5){
      rect(utr5[["start"]],1.5, utr5[["end"]],2.5,
           col="orange",
           border="black")
    })
  }
  
  
  if(nrow(intergenic_coordinates)>0){
    apply(intergenic_coordinates,1, function(inter){
      rect(inter[["start"]],1.5, inter[["end"]],2.5,
           col="black",
           border="black")
    })
  }
  
  if(nrow(DMRs_of_interest)>0){
    coords <- apply(DMRs_of_interest,1,function(dmr){
      if(dmr[["is.hyper"]]=="hyper"){col="red"}
      else{col="blue"}
      rect(dmr[["start"]],1.1,dmr[["end"]],1.3,
           col=col,
           border=col)
      coord <- as.numeric(c(dmr[["start"]],dmr[["end"]]))
      
    })
    colnames(coords)<-DMRs_of_interest[["DMR_id"]]
    coords <- apply(coords,2,mean)
    text(coords,rep(1.2,length(coords)),names(coords),cex = 0.5)
  }
  
  
  
  plot_coordinates <<- list(introns_coordinates = introns_coordinates,
                            exons_coordinates = exons_coordinates,
                            utr3_coordinates = utr3_coordinates,
                            utr5_coordinates = utr5_coordinates,
                            p2000_coordinates = p2000_coordinates,
                            intergenic_coordinates = intergenic_coordinates,
                            probes_coordinates = sub_epic,
                            DMRs_of_interest = DMRs_of_interest)
  
  
  
  
  
}




