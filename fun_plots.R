##################### Plot probes position along a selected gene given a features object #############

plot_selected_probes = function(selected_gene = gene, expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = probes_index){
  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))] #### On ne garde en données que les individus pour lesquels on a méthylation ET 
  
  selected_expr_values = sort(expr_data[selected_gene,]) ##### On sélectionne les données d'expression du gène sélectionné
  selected_probes = intersect(features$probes[[selected_gene]],rownames(meth_data)) #### On sélectionne les probes se trouvants sur le gène sélectionné
  tmp = meth_lusc$platform[selected_probes,]
  
  plot(tmp$Start,rep(1,length(tmp$Start)),yaxt = "n",ylab ="",xlab = "Probes physical coordinates along the gene (in bp)",main=paste("Probes repartition along",selected_gene,sep=" "),pch=16)
  text(tmp$Start,rep(1.1,length(tmp$Start)), labels=selected_probes,cex = 0.8,srt = 80)
  text(features$features[selected_gene,"TSS"],0.95, labels="TSS",cex = 0.9)
  points(features$features[selected_gene,"TSS"],1,col="red",pch = 4)
  
  
  
}


######################## Plotting transcription/expression level & methylation for a selected gene among pre-created features############################



plot_selected_gene  = function(selected_gene = "DKKL1", expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = probes_index){
  
  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))] #### On ne garde en données que les individus pour lesquels on a méthylation ET 
  
  selected_expr_values = sort(expr_data[selected_gene,]) ##### On sélectionne les données d'expression du gène sélectionné
  selected_probes = probes_index[[selected_gene]] #### On sélectionne les probes se trouvants sur le gène sélectionné
  methvals_of_interest = meth_data[intersect(selected_probes,rownames(meth_data)),order(match(colnames(meth_data),names(selected_expr_values)))]
  ## ci-dessus on prend les données de méthylations émises par nos sondes d'intérêts, dans une df triée selon le niveau d'expression du gène sélectionné par patient.
  
  
  
  
  layout(matrix(1:2,1), respect=TRUE)
  expression_plot = plot(selected_expr_values, xlab = paste("Patient index ( ordered by level of expression)"), ylab = paste("Expression level"), main = "Expression/Transcription")
  
  heatmap = image(methvals_of_interest ,ylab = "Patient index ( ordered by level of expression)", main = "Methylation", Rowv = NA, Colv = NA,axes=FALSE)
  #axis(side = 1, at=seq(0,1,1/(nrow(methvals_of_interest)-1)), labels = rownames(methvals_of_interest), tick = FALSE, las=2)
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
}
