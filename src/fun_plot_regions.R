
plot_selected_regions<- function(selected_gene = "ALDH3B1", features_list = subfeatures, pf = platform, EPIMEDepic = epic, probes_index = feat_indexed_probes){
  
  ## get data
  
    pf_gene <- platform[which(platform[,"gene"]==selected_gene),2:4]
    tmp_probes <- feat_indexed_probes[[selected_gene]]
    epic_tmp <- epic[tmp_probes,c("start","end")]
    sub_epic <- epic_tmp[intersect(rownames(epic_tmp), rownames(meth_lusc$data)),]
    TSS <- features_list[selected_gene,"TSS"]
    
    ## get subset of coordinates
    
    introns_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTRON"),]
    exons_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="CDS"),]
    utr3_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="3UTR"),]
    utr5_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="5UTR"),]
    intergenic_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTER"),]
    p2000_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="P2000"),]
    
    ## plot probes + TSS
    
    plot(sub_epic$start, rep(3,length(sub_epic$start)), pch=3, xlim=c(TSS-10000,TSS+10000), cex=1, yaxt="n",
         main = paste0(selected_gene,": Probes repartition among regions (nprobes = ",nrow(sub_epic),")"),
         xlab="Coordinates (in bp)",
         ylab="",
         ylim=c(0,5))
    
    legend("left", c("P2000", "Exons", "Introns","3'UTR","5'UTR","Inter"),inset=c(-0.12,0), xpd = TRUE, pch=15, col=c("grey","blue","red","green","orange","black"),bty="n")
    
    text(sub_epic$start,rep(4,length(sub_epic$start)), labels = rownames(sub_epic) ,cex = 0.8,srt = 80)
    text(TSS,0.2,labels = paste0("TSS : ", TSS),cex = 0.7)
    
    points(TSS, 0.5, pch=9, col="red")
    abline(h=0.5, lty=3, lwd = 1)
    
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
               col="blue",
               border="black")
      })
    }
    
    
    if(nrow(introns_coordinates)>0){
      apply(introns_coordinates,1, function(intron){
        rect(intron[["start"]],1.5, intron[["end"]],2.5,
               col="red",
               border="black")
      })
    }
    
    if(nrow(utr3_coordinates)>0){
      apply(utr3_coordinates,1, function(utr3){
        rect(utr3[["start"]],1.5, utr3[["end"]],2.5,
             col="green",
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
    
    
    ret = list(introns_coordinates = introns_coordinates,
               exons_coordinates = exons_coordinates,
               utr3_coordinates = utr3_coordinates,
               utr5_coordinates = utr5_coordinates,
               p2000_coordinates = p2000_coordinates,
               intergenic_coordinates = intergenic_coordinates,
               probes_coordinates = sub_epic)
    
    return(ret)
    
  
  
}




