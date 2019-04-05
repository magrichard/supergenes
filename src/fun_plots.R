##################### Plot probes position along a selected gene given a features object #############

plot_selected_probes <- function(selected_gene = "ALDH3B1",binlist=c("bin1","bin2","bin3","bin4","bin5","bin6"), meth_data = meth_lusc$data, trscr_study = trscr_lusc){

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
  
  
  
  
  tmp = means_per_bins_per_genes_per_patient[startsWith(rownames(means_per_bins_per_genes_per_patient),selected_gene),]
  xlabs = rownames(tmp)
  boxplot(t(tmp),xaxt="n",main = "Distribution of mean Mval per bin")
  axis(side = 1, at=1:6, labels = xlabs, tick = FALSE, las=2)
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
}


######################## Plotting transcription/expression level & methylation for a selected gene among pre-created features############################



plot_selected_gene  = function(selected_gene = selected_gene, expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = feat_indexed_probes, ...){
  
  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))]  
  selected_expr_values = sort(expr_data[selected_gene,]) 
  selected_probes = feat_indexed_probes[[selected_gene]] 
  methvals_of_interest = meth_data[intersect(selected_probes,rownames(meth_data)),order(match(colnames(meth_data),names(selected_expr_values)))]
  
  
  
  
  
  layout(matrix(1:2,1), respect=TRUE)
  tmp = cbind(selected_expr_values,1:length(selected_expr_values))
  expression_plot = plot(tmp, ylab = paste("Patient index ( ordered by level of expression)"), xlab = paste("Expression level"), main = "Expression/Transcription",xlim=c(0,20))
  
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(100)
  
  heatmap = image(methvals_of_interest, Rowv = NA, Colv = NA, axes=FALSE, col=cols, main="Methylation")
  axis(side = 1, at=seq(0,1,1/(nrow(methvals_of_interest)-1)), labels = rownames(methvals_of_interest), tick = FALSE, las=2)
  legend("left", c("hypo", "neutral", "hyper"), xpd = TRUE, pch=15, inset = c(-0.35,-0.25), col=c("green","black","red"),bty="n")
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
}





