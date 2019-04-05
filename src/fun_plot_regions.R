plot_selected_regions(selected_gene = "ALDH3B1", features_list = subfeatures, platform = platform, EPIMEDepic = epic, probes_index = feat_indexed_probes){
  
  
  
}

  pf_gene <- platform[which(platform[,"gene"]==selected_gene),2:4]
  tmp_probes <- feat_indexed_probes[[selected_gene]]
  epic_tmp <- epic[tmp_probes,c("start","end")]
  sub_epic <- epic_tmp[intersect(rownames(epic_tmp), rownames(meth_lusc$data)),]
  TSS <- features_list[selected_gene,"TSS"]
  
  introns_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTRON"),]
  exons_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="CDS"),]
  utr3_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="3UTR"),]
  utr5_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="5UTR"),]
  intergenic_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTER"),]
  p2000_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="P2000"),]

plot(sub_epic$start, rep(1,length(sub_epic$start)),pch=16,xlim=c(TSS-10000,TSS+10000),cex=0.75)
points(TSS,0.92,pch=9,col="red")
abline(h=0.92, lty=3, lwd = 1)

    if(nrow(p2000_coordinates) > 0) {
  points(p2000_coordinates[,"start"],rep(0.95,nrow(p2000_coordinates)),pch=0)
    }
  
  if(nrow(exons_coordinates)>0){
    points(exons_coordinates[,"start"],rep(0.93,nrow(exons_coordinates)),pch=3,col="blue")
  }
  
  if(nrow(introns_coordinates)>0){
    points(introns_coordinates[,"start"],rep(0.95,nrow(introns_coordinates)),pch=3,col="green")
  }
  if(nrow(utr3_coordinates)>0){
    points(utr3_coordinates[,"start"],rep(0.95,nrow(utr3_coordinates)))
  }
  if(nrow(utr5_coordinates)>0){
    points(utr5_coordinates[,"start"],rep(0.95,nrow(intergenic_coordinates)))
  }
  if(nrow(intergenic_coordinates)>0){
    points(intergenic_coordinates[,"start"],rep(0.95,nrow(intergenic_coordinates)),pch=0)
  }



plot(sub_epic$start,rep(1,length(sub_epic$start)),yaxt="n",xaxt="n",xlim=c(TSS-2500,TSS+2500),ylab ="", xlab = "", main="Probes along the gene (in bp)",pch=16)
text(sub_epic$start,rep(1.1,length(sub_epic$start)), labels=selected_probes,cex = 0.8,srt = 80)
bins_coord = c(TSS-2500,TSS-1000,TSS-500,TSS,TSS+500,TSS+1000,TSS+2500)
sapply(bins_coord, function(bin){abline(v=bin,lty=3)})
axis(side = 1, at=bins_coord, labels = c("-2500","-1000","-500","TSS","+500","+1000","+2500"), tick = FALSE, las=2)
mtext(paste(noquote(tmp)), cex=1, line=-3)

plot(sub_epic$start, rep(1,length(sub_epic$start)))
