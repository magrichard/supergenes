
################### plot mean vs sd vals per probe ######################



#features <- get_features(penda_superup_deregulated, study = trscr_lusc, up_str = 7500, dwn_str = 7500)
#layout(matrix(1:2,1), respect=TRUE)

plot_mean_vs_sd <- function(genes_list = feat_indexed_probes,
                            data = meth_tumoral,...){
  
  
  
  poi <-  unique(unlist(genes_list))
  meth <- data[intersect(rownames(data),poi),]
  
  rsd_probes <- epimedtools::monitored_apply(meth,1, rsd,na.rm=T)
  sd_probes <- epimedtools::monitored_apply(meth,1, sd,na.rm=T)
  means_probes <- epimedtools::monitored_apply(meth,1, mean,na.rm=T)
  
  tmp_stats <- cbind(means_probes[!is.na(means_probes)],sd_probes[!is.na(sd_probes)],rsd_probes[!is.na(rsd_probes)])
  stats <- tmp_stats[order(tmp_stats[,1]),]
  colnames(stats)<- c("means","sd","rsd")
  
  
  plot(stats[,"means"],stats[,"sd"],
       xlab = paste0("mean per probe \n (nprobes = ",nrow(stats)," + ",nrow(meth)-nrow(stats)," NAs )"),
       ylab="sd per probe",...)
  
  return(stats)
}
















##################### Plot probes position along a selected gene given a features object #############

plot_genes_probes <- function(selected_gene = "ALDH3B1",binlist=c("bin1","bin2","bin3","bin4","bin5","bin6"), meth_data = meth_lusc$data, trscr_study = trscr_lusc, distributions_per_bins = means_per_bins_per_genes_per_patient){

  if(trscr_lusc$platform[selected_gene,"strand"]=="+") {tmp = "3'->5'"} else {tmp = "3'->5'"}
  
  
  
  selected_probes = intersect(feat_indexed_probes[[selected_gene]],rownames(meth_data))
  sub_epic<-feat_indexed_probes_bin[[selected_gene]][[1]]
  TSS <- features[selected_gene,"TSS"]
  
  
  layout(matrix(1:2,1), respect=TRUE)
  
  plot(sub_epic$start,rep(1,length(sub_epic$start)),yaxt="n",xaxt="n",xlim=c(TSS-2500,TSS+2500),ylab ="", xlab = "", main="Probes along the gene (in bp)",pch=16)
  text(sub_epic$start,rep(1.1,length(sub_epic$start)), labels=selected_probes,cex = 0.8,srt = 80)
  bins_coord = c(TSS-2500,TSS-1000,TSS-500,TSS,TSS+500,TSS+1000,TSS+2500)
  sapply(bins_coord, function(bin){abline(v=bin,lty=3)})
  axis(side = 1, at=bins_coord, labels = c("-2500","-1000","-500","TSS","+500","+1000","+2500"), tick = FALSE, las=2)
  mtext(paste(noquote(tmp)), cex=1, line=-3)
  
  
  
  
  tmp = distributions_per_bins[startsWith(rownames(distributions_per_bins),selected_gene),]
  xlabs = rownames(tmp)
  boxplot(t(tmp),xaxt="n",main = "Distribution of mean Mval per bin")
  axis(side = 1, at=1:6, labels = xlabs, tick = FALSE, las=2)
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
}


######################## Plotting transcription/expression level & methylation for a selected gene among pre-created features############################



plot_gene_meth  = function(selected_gene = selected_gene, expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = feat_indexed_probes, ...){
  
  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))]  
  selected_expr_values = sort(expr_data[selected_gene,]) 
  selected_probes = feat_indexed_probes[[selected_gene]] 
  methvals_of_interest = meth_data[intersect(selected_probes,rownames(meth_data)),order(match(colnames(meth_data),names(selected_expr_values)))]
  
  
  
  
  
  layout(matrix(1:2,1), respect=TRUE)
  tmp = cbind(selected_expr_values,1:length(selected_expr_values))
  expression_plot = plot(tmp, ylab = paste("Patient index ( ordered by level of expression)"), xlab = paste("Expression level"), main = "Expression/Transcription",xlim=c(0,20))
  
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(100)
  
  breaks <- c(seq(0,0.33,length=35),seq(0.34,0.66,length=33),seq(0.67,1,length=33))
  
  heatmap = image(methvals_of_interest, axes=FALSE, col=cols, main="Methylation",breaks = breaks)
  axis(side = 1, at=seq(0,1,1/(nrow(methvals_of_interest)-1)), labels = rownames(methvals_of_interest), tick = FALSE, las=2)
  legend("left", c("hypo", "neutral", "hyper"), xpd = TRUE, pch=15, inset = c(-0.35,-0.25), col=c("green","black","red"),bty="n")
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
  

}


############Plot genes with regions##############



  plot_genes_regions<- function(selected_gene = "ALDH3B1",
    features_list = features,
    DMRtable = DMR_table,
    cgi_coordinates = cgi,
    DMR_map, pf = platform,
    epic_850k = epic,
    epic_27k = epic27k,
    epic_450k = epic450k,
    probes_index = feat_indexed_probes,
    window=c(10000,10000)){
  
    ## get nearest DMR
    
    DMRs_of_interest<-DMRtable[which(DMRtable[,"closest_gene_name"]==selected_gene),]
    #DMRs_of_interest<- DMR_table[intersect(DMR_table[,"DMR_id"],names(DMR_map[[selected_gene]])),]
    
    
    
    ## get data
    upstream = window[2]
    downstream = window[1]
    pf_gene <- platform[which(platform[,"gene"]==selected_gene),2:4]
    tmp_probes <- feat_indexed_probes[[selected_gene]]
    epic_tmp <- epic_450k[tmp_probes,c("start","End")]
    sub_epic <- epic_tmp[intersect(rownames(epic_tmp), rownames(meth_lusc$data)),]
    sub_epic <- cbind(sub_epic,cgi_coordinates[intersect(rownames(sub_epic),rownames(cgi_coordinates)),])
    TSS <- features_list[selected_gene,"TSS"]
    
    ## get subset of coordinates
    
    introns_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTRON"),]
    exons_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="CDS"),]
    utr3_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="3UTR"),]
    utr5_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="5UTR"),]
    intergenic_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="INTER"),]
    p2000_coordinates <- pf_gene[which(pf_gene[,"genomic_feature"]=="P2000"),]
  
  ## get colors
  
    cols<-sapply(rownames(sub_epic), function(probe){
    
      if(is.element(probe,rownames(epic450k))){
        col ="darkcyan"
      }
      
      if(is.element(probe,rownames(epic27k))){
        
        col = "chartreuse4"
      }
      if(is.element(probe,rownames(epic27k)) & is.element(probe,rownames(epic450k))){
        
        col = "red"
      }
      if(!is.element(probe,rownames(epic27k)) & !is.element(probe,rownames(epic450k))){
        
        col = "black"
      }
      
      
    return(col)
  })
  
  ## plot probes + TSS
  
  plot(sub_epic$start, rep(3,length(sub_epic$start)), pch=19, xlim=c(TSS-downstream,TSS+upstream), cex=1.2, yaxt="n",
       main = paste0(selected_gene,": Probes repartition among regions (nprobes = ",nrow(sub_epic),")"),
       xlab="Coordinates (in bp)",
       ylab="",
       ylim=c(0,5),
       col = cols)
  
  legend("left", c("P2000", "Exons", "Introns","3'UTR","5'UTR","Inter"),inset=c(-0.12,0), xpd = TRUE, pch=15, col=c("grey","blue","red","green","orange","black"),bty="n",title = "regions")
  legend("topright", c("450k","27k","both","850k only"),inset=c(-0.05,0), xpd=TRUE, pch = 19, col = c("darkcyan","chartreuse4","red","black"), bty="n",title ="Chip")
  #text(sub_epic$start,rep(4,length(sub_epic$start)), labels = rownames(sub_epic) ,cex = 0.8,srt = 80) ##probes label
  text(TSS,0.2,labels = paste0("TSS : ", TSS),cex = 0.7)
  
  points(TSS, 0.5, pch=9, col="red")
  abline(h=0.5, lty=3, lwd = 1)
  
  
  
  ##plot CGIs
  
  for (i in 1:21){
    
    rect(as.numeric(as.character(sub_epic[i,"cgi_start"])),2.7, as.numeric(as.character(sub_epic[i,"cgi_end"])),2.6,
         density = 5)
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


  #################heatmap######################
  
  meth_heatmap <-function(data = means, dendrogram = "none", Rowv = NULL, Colv = NULL,cols = NULL, ...){
    
    if(!exists("main")){
      main = deparse(substitute(data))
    }
    
    main=deparse(substitute(data))
    data = data[,-7]
    
    if(is.null(cols)){
      colors=c("green", "black", "red")
      cols = colorRampPalette(colors)(100)
    }
    
    foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, trace="none", col=cols, mar=c(10,5), useRaster=TRUE ,...)
    
    
}
  
