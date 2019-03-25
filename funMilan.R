catalog_selected_genes =  function(selected_features = features){
  
  available_genes = list(catalog="catalog", n_genes = "n_genes")
  
  available_genes[["catalog"]] = noquote(rownames(features[!is.na(features$meth_cluster),]))
  available_genes [["n_genes"]] = noquote(paste("Genes selected :", length(available_genes[[1]]),sep = " "))
  warning("Note that quotes were removed for clarity, please select your genes using it for the plotting function.")
  return(available_genes)
  
}




trscr_lusc = readRDS("~/projects/supergenes/data/tcga_studies/study_TCGA-LUSC_trscr.rds")



plot_selected_gene  = function(selected_gene = "DKKL1", expr_data = expr_trscr_lusc$data, meth_data = study_dmeth_tcga_lusc, probes_index = feat_indexed_probes){
  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))] #### On ne garde en données que les individus pour lesquels on a méthylation ET 
  
  selected_expr_values = sort(expr_data[selected_gene,]) ##### On sélectionne les données d'expression du gène sélectionné
  selected_probes = feat_indexed_probes[[selected_gene]] #### On sélectionne les probes se trouvants sur le gène sélectionné
  methvals_of_interest = study_dmeth_tcga_lusc[intersect(selected_probes,rownames(study_dmeth_tcga_lusc)),order(match(colnames(study_dmeth_tcga_lusc),names(selected_expr_values)))]
  ## ci-dessus on prend les données de méthylations émises par nos sondes d'intérêts, dans une df triée selon le niveau d'expression du gène sélectionné par patient.
  
  
  
  
  layout(matrix(1:2,1), respect=TRUE)
  expression_plot = plot(selected_expr_values, xlab = paste("Patient index ( ordered by level of expression)"), ylab = paste("Expression level"), main = "Expression/Transcription")
  
  heatmap = image(methvals_of_interest ,ylab = "Patient index ( ordered by level of expression)", main = "Methylation", Rowv = NA, Colv = NA,axes=FALSE)
  axis(side = 1, at=seq(0,1,1/(nrow(methvals_of_interest)-1)), labels = rownames(methvals_of_interest), tick = FALSE, las=2)
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
}



