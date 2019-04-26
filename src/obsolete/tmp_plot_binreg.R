



plot_binreg<- function(selected_gene = "ANLN",
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
  
  strand <- features_list[selected_gene,"strand"]
  TSS <- features_list[selected_gene,"TSS"]
  upstream = window[2]
  downstream = window[1]
  
  pf_gene <- pf[which(pf[,"gene"]==selected_gene),2:4]
  tmp_probes <- feat_indexed_probes[[selected_gene]]
  sub_epic_tmp <- epimed[tmp_probes,c("start","end")]
  sub_epic_tmp <- sub_epic_tmp[intersect(rownames(sub_epic_tmp), rownames(meth_lusc$data)),]
  sub_epic <- cbind(sub_epic_tmp,cgi_table[intersect(rownames(sub_epic_tmp),rownames(cgi_table)),])
  
  chip_index_of_interest <- chips_index[intersect(rownames(sub_epic),rownames(chips_index)),]
  
  
  ## get subset of coordinates
  if(strand == "+"){ limit <- "start"} else (limit <- "end")
    
  
  introns_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTRON"),]
  exons_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="CDS"),]
  utr3_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="3UTR"),]
  utr5_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="5UTR"),]
  intergenic_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTER"),]
  
  ## keep only regions downstream for strand = "+"
  
  if(strand == "+"){
  introns_used <- introns_coordinates[which(introns_coordinates[,"end"] > TSS),]
  exons_used <-exons_coordinates[which(exons_coordinates[,"end"] > TSS),]
  utr3_used <-utr3_coordinates[which(utr3_coordinates[,"end"] > TSS),]
  utr5_used <-utr5_coordinates[which(utr5_coordinates[,"end"] > TSS),]
  intergenic_used <-intergenic_coordinates[which(intergenic_coordinates[,"end"] > TSS),]
  
  limits <- c(TSS,TSS-500,TSS-2500,TSS-1000,TSS-500)
  }
  
  ## keep only regions upstream for strand = "-"
  
  if(strand == "-"){
    introns_used <- introns_coordinates[which(introns_coordinates[,"end"] < TSS),]
    exons_used <-exons_coordinates[which(exons_coordinates[,"end"] < TSS),]
    utr3_used <-utr3_coordinates[which(utr3_coordinates[,"end"] < TSS),]
    utr5_used <-utr5_coordinates[which(utr5_coordinates[,"end"] < TSS),]
    intergenic_used <-intergenic_coordinates[which(intergenic_coordinates[,"end"] < TSS),]
    
    limits <- c(TSS,TSS+500,TSS+2500,TSS+1000,TSS+500)
  }
  
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
  
  plot(sub_epic$start, rep(3,length(sub_epic$start)), pch=19, xlim=c(TSS-downstream,TSS+upstream), cex=1.1, yaxt="n",
       main = paste0(selected_gene,": Probes repartition among regions & bins \n (nprobes = ",nrow(sub_epic),")"),
       xlab="Coordinates (in bp)",
       ylab="",
       ylim=c(0,5),
       col = cols)
  
  legend("left", c("Exons", "Introns","3'UTR","5'UTR","Inter","hypoDMR","hyperDMR","CGI"),
         inset=c(-0.12,0),
         xpd = TRUE,
         pch=15,
         col=c("darkslategray1","firebrick","darkseagreen","orange","black","blue","red","gray1"),
         title = "Region")
  
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
  abline(h=0.5, lty=1, lwd = 1)
  
  
  sapply(limits, function(limit){abline(v=limit,lty=3)})
  
  
  
  ##plot CGIs
  
  for (i in 1:nrow(sub_epic)){
    
    rect(as.numeric(as.character(sub_epic[i,"cgi_start"])),2.7, as.numeric(as.character(sub_epic[i,"cgi_end"])),2.6,
         col ="gray1")
  }
  
  
  
  
  ## plot regions
  
  
  
  
  
  if (nrow(exons_used)>0){
    apply(exons_used,1, function(exon){
      rect(exon[["start"]],1.5, exon[["end"]], 2.5,
           col="darkslategray1",
           border="black")
    })
  }
  
  
  if(nrow(introns_used)>0){
    apply(introns_used,1, function(intron){
      rect(intron[["start"]],1.5, intron[["end"]],2.5,
           col="firebrick",
           border="black")
    })
  }
  
  if(nrow(utr3_used)>0){
    apply(utr3_used,1, function(utr3){
      rect(utr3[["start"]],1.5, utr3[["end"]],2.5,
           col="darkseagreen",
           border="black")
      
    })
  }
  
  
  
  if(nrow(utr5_used)>0){
    apply(utr5_used,1, function(utr5){
      rect(utr5[["start"]],1.5, utr5[["end"]],2.5,
           col="orange",
           border="black")
    })
  }
  
  
  if(nrow(intergenic_used)>0){
    apply(intergenic_used,1, function(inter){
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
  
  
  
  plot_coordinates <<- list(introns_used = introns_used,
                            exons_used = exons_used,
                            utr3_used = utr3_used,
                            utr5_used = utr5_used,
                            intergenic_used = intergenic_used,
                            probes_coordinates = sub_epic,
                            DMRs_of_interest = DMRs_of_interest)
  
  
  
  
  
}
