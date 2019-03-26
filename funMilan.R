################################### Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################


get_features <-function(targeted_genes,index_pathology = 2 ,study){   #### In my uses-case, 1 for LUAD & 2 for LUSC 
  print("Creating a list of features...")
  targeted_genes = targeted_genes[targeted_genes[,index_pathology]==1,]
  features = study$platform[intersect(rownames(study$platform),rownames(targeted_genes)),]
  return(features)
  
}


### Ex : get_features(penda_superconserved,index_pathology = 2,trscr_lusc)

############################ Visualisation tool for features created by get_features() function##########################################


catalog_selected_genes =  function(selected_features = sub_features){
  
  available_genes = list(catalog="catalog", n_genes = "n_genes")
  
  available_genes[["catalog"]] = noquote(rownames(selected_features))
  available_genes [["n_genes"]] = noquote(paste("Genes selected :", length(available_genes[[1]]),sep = " "))
  warning("Note that quotes were removed for clarity, please select your genes using it for the plotting function.")
  return(available_genes)
  
}

##################### Création des bins #############

desc_selected_gene = function(selected_gene = gene, expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = feat_indexed_probes){

  tmp = probes_index[[selected_gene]]
  tmp2 = meth_data[intersect(tmp,meth_data[]),]
  


  

  

  

  
  


  

}


######################## Plotting transcription/expression level & methylation for a selected gene among pre-created features############################



plot_selected_gene  = function(selected_gene = "DKKL1", expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = feat_indexed_probes){

  
  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))] #### On ne garde en données que les individus pour lesquels on a méthylation ET 
  
  selected_expr_values = sort(expr_data[selected_gene,]) ##### On sélectionne les données d'expression du gène sélectionné
  selected_probes = feat_indexed_probes[[selected_gene]] #### On sélectionne les probes se trouvants sur le gène sélectionné
  methvals_of_interest = meth_data[intersect(selected_probes,rownames(meth_data)),order(match(colnames(meth_data),names(selected_expr_values)))]
  ## ci-dessus on prend les données de méthylations émises par nos sondes d'intérêts, dans une df triée selon le niveau d'expression du gène sélectionné par patient.
  
  
  
  
  layout(matrix(1:2,1), respect=TRUE)
  expression_plot = plot(selected_expr_values, xlab = paste("Patient index ( ordered by level of expression)"), ylab = paste("Expression level"), main = "Expression/Transcription")
  
  heatmap = image(methvals_of_interest ,ylab = "Patient index ( ordered by level of expression)", main = "Methylation", Rowv = NA, Colv = NA,axes=FALSE)
  axis(side = 1, at=seq(0,1,1/(nrow(methvals_of_interest)-1)), labels = rownames(methvals_of_interest), tick = FALSE, las=2)
  
  mtext(paste(noquote(selected_gene),"Analysis",sep= " "), outer=TRUE,  cex=2, line=-2)
  
}



