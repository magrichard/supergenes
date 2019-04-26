
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

plot_genes_probes <- function(selected_gene = "ALDH3B1",
                              binlist=c("bin1","bin2","bin3","bin4","bin5","bin6"),
                              meth_data = meth_lusc$data,
                              trscr_study = trscr_lusc,
                              distributions_per_bins = means_per_bins_per_genes_per_patient){

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



plot_gene_meth  = function(selected_gene = selected_gene,
                           expr_data = trscr_lusc$data,
                           meth_data = meth_lusc$data,
                           probes_index = feat_indexed_probes,
                           ...){
  
  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))]  
  selected_expr_values = sort(expr_data[selected_gene,]) 
  selected_probes = feat_indexed_probes[[selected_gene]] 
  methvals_of_interest = meth_data[intersect(selected_probes,rownames(meth_data)),order(match(colnames(meth_data),names(selected_expr_values)))]
  
  
  
  
  
  layout(matrix(1:2,1), respect=TRUE)
  tmp = cbind(selected_expr_values,1:length(selected_expr_values))
  expression_plot = plot(tmp, ylab = paste("Patient index ( ordered by level of expression)"), xlab = paste("Expression level"), main = "Expression/Transcription",xlim=c(0,20))
  
  colors=c("cyan", "black", "red")
  cols = colorRampPalette(colors)(100)
  
  breaks <- c(seq(0,0.33,length=35),seq(0.34,0.66,length=33),seq(0.67,1,length=33))
  
  heatmap = image(methvals_of_interest, axes=FALSE, col=cols, main="Methylation",breaks = breaks)
  axis(side = 1, at=seq(0,1,1/(nrow(methvals_of_interest)-1)), labels = rownames(methvals_of_interest), tick = FALSE, las=2)
  legend("left", c("hypo", "neutral", "hyper"), xpd = TRUE, pch=15, inset = c(-0.35,-0.25), col=c("cyan","black","red"),bty="n")
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
  

}


############Plot genes with regions##############



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


  #################heatmap######################
  
  meth_heatmap <-function(data = means, dendrogram = "none", Rowv = NULL, Colv = NULL,cols = NULL, ...){
    
    if(!exists("main")){
      main = deparse(substitute(data))
    }
    
    main=deparse(substitute(data))
    data = data[,-which(colnames(data)=="overall")]
    
    if(is.null(cols)){
      colors=c("cyan", "black", "red")
      cols = colorRampPalette(colors)(100)
    }
    
    foo = gplots::heatmap.2(data,
                            Rowv=Rowv,
                            Colv=Colv,
                            dendrogram=dendrogram,
                            trace="none",
                            col=cols,
                            mar=c(10,5),
                            useRaster=TRUE,
                            ...)
    
    
}
  
############## plot gene expr along with its associated DMRs methylation ############
  
  
  
  plot_gene_dmr<- function(selected_gene = "OTX1",
                           expr_data = trscr_lusc$data,
                           meth_study = meth_lusc,
                           dmrs_table = DMR_table,
                           DMR_map  = DMRs100k,
                           fun = mean){
    
    
    
    tmp_DMR <- dmrs_table[names(DMR_map[[selected_gene]]),c("chr","start","end","is.hyper")]
    DMRs_of_interest<- tmp_DMR[order(tmp_DMR[,"start"]),]
    
    if(nrow(DMRs_of_interest) == 0) {stop(paste0("There is no DMRs indexed for ",selected_gene,". Check your DMR_map object"))}
    
    probes_on_chr <- epic[which(as.character(epic[,1]) == DMRs_of_interest[["chr"]][1]),]
    
    
    expr_values <- sort(expr_data[selected_gene,intersect(colnames(expr_data),colnames(meth_study$data))])
    meth_values <- meth_study$data[,intersect(colnames(expr_data),colnames(meth_study$data))]
    
    
    
    
    
    
    
    
    
    tmp_DMR <- dmrs_table[names(DMR_map[[selected_gene]]),c("chr","start","end","is.hyper")]
    DMRs_of_interest<- tmp_DMR[order(tmp_DMR[,"start"]),]
    probes_on_chr <- epic[which(as.character(epic[,1]) == DMRs_of_interest[["chr"]][1]),]
    
    print("Indexing methylation values per DMR...")
    vals_per_dmr<-epimedtools::monitored_apply(t(t(rownames(DMRs_of_interest))),1, function(i){
      
      DMR<-DMRs_of_interest[i,]   
      
      
      selected_probes <- probes_on_chr[which(probes_on_chr[,"start"] >= DMR[["start"]] & probes_on_chr[,"end"] <= DMR[["end"]]),]
      selected_meth_values <- meth_values[intersect(rownames(selected_probes),rownames(meth_values)),]
      
      
      if(class(selected_meth_values) != "numeric"){
        DMR_vals <- apply(selected_meth_values,2,fun,na.rm=T)}         #### Check if there is more than 2 selected probes
      
      else{DMR_vals <- rep(NA,length(selected_meth_values))}
      
      
      
      
      #### Else, returning raw values if the function is mean, else return a vector of NAs
      
      
      return(DMR_vals)
      
      
    })
    
    rownames(vals_per_dmr)<- colnames(meth_values)
    sorted_vals = vals_per_dmr[order(match(rownames(vals_per_dmr),names(expr_values))),]
    colnames(sorted_vals)<- rownames(DMRs_of_interest)
    
    
    cols1 <- epimedtools::monitored_apply(t(t(names(expr_values))),1, function(i){
      if(meth_study$exp_grp[i,14] == "tumoral"){col<- "red"}
      else {col<- "blue"}})
    
    print("Plotting...")
    
    layout(matrix(1:2,1), respect=TRUE)
    
    
    tmp = cbind(expr_values,1:length(expr_values))
    expression_plot = plot(tmp,
                           ylab = paste("Patient index ( ordered by level of expression)"),
                           xlab = paste("Expression level"),
                           main = "Expression/Transcription",
                           col = cols1,
                           xlim=c(0,20))
    legend("topleft",
           c("tumoral","normal"),
           xpd = TRUE,
           pch=20,
           col=c("red","blue"))
    
    
    colors=c("cyan", "black", "red")
    cols = colorRampPalette(colors)(100)
    breaks <- c(seq(0,0.33,length=35),seq(0.34,0.66,length=33),seq(0.67,1,length=33))
    
    image(t(sorted_vals),
          col=cols,
          breaks = breaks,
          axes=FALSE,
          main = paste0("# of DMRs : ",ncol(vals_per_dmr)))
    
    mtext(paste(noquote(selected_gene),"DMRs vizualisation",sep= " "), outer=TRUE,  cex=2, line=-2)
    
    
    return(sorted_vals)
    
    
    
}
  
  #table(unlist(lapply(DMRs100k, length)))  
  
######## Plot bin & regions ###############
  
  
  
  
  
  
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
      introns_used <- introns_coordinates[which(introns_coordinates[,"start"] < TSS),]
      exons_used <-exons_coordinates[which(exons_coordinates[,"start"] < TSS),]
      utr3_used <-utr3_coordinates[which(utr3_coordinates[,"start"] < TSS),]
      utr5_used <-utr5_coordinates[which(utr5_coordinates[,"start"] < TSS),]
      intergenic_used <-intergenic_coordinates[which(intergenic_coordinates[,"start"] < TSS),]
      
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
  
  
######## boxplot results ( compare 2 outputs of subset_vals_per_bins function) ##########
  
boxplot_res <- function(values_h = means_h,
                        values_t = means_t,
                        main_h = "healty tissues",
                        main_t="tumoral tissues",
                        levels = TRUE){
    
    layout(matrix(2:1,1), respect=TRUE)
    
    boxplot(values_h,main = main_h)
    means_h<-apply(values_h,2,mean,na.rm=T)
    points(means_h,col="red",pch = 19,cex = 0.5)
    
    boxplot(values_t,main = main_t)
    means_t<-apply(values_t,2,mean,na.rm=T)
    points(means_t,col="red",pch = 19,cex= 0.5)
    
    
    
    mtext(paste0(deparse(substitute(values_h))," vs ",deparse(substitute(values_t))), outer=TRUE,  cex=2, line=-2)
    
    if(isTRUE(levels)){
    print(noquote(c("Factors for graphic h :",paste0(colnames(values_h)))))
    print(noquote(c("Factors for graphic t :",paste0(colnames(values_t)))))
    }
}
    
###########  sd vs mean ###########

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
       xlab = paste0("mean per probe (nprobes = ",nrow(stats)," + ",nrow(meth)-nrow(stats)," NAs )"),
       ylab="sd per probe",...)
  
  return(stats)
}
  

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
  
  plot(sub_epic$Start, rep(1,length(sub_epic$Start)), pch=19, xlim=c(closest_cgi[["center"]]-window[1],closest_cgi[["center"]]+window[2]), cex=1.1, yaxt="n",
       main = paste0("nprobes = ",nrow(sub_epic),", ncgi = ",nrow(cgi_of_interest)),
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
  
  breaks <- c(seq(0,0.33,length=35),
              seq(0.34,0.66,length=33),
              seq(0.67,1,length=33))
  
  image(methvals_of_interest,
        axes=FALSE,
        col=cols,
        main=paste0("nprobes indexed : ",nrow(methvals_of_interest)),
        breaks = breaks,
        ylab="Patients",
        xlab="probes")
  
  mtext(paste(noquote(selected_gene),"Methylation profile",sep= " "), outer=TRUE,  cex=2, line=-2)
  
  
} 
  
  