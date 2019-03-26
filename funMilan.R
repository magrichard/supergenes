################################### Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################


get_features <- function(targeted_genes, index_pathology = 2, study, up_str = 2500, dwn_str = 2500, nb_probe_min = 1){   #### In my uses-case, 1 for LUAD & 2 for LUSC 
  print("Creating a list of features...")
  targeted_genes = targeted_genes[targeted_genes[,index_pathology]==1,]
  features = study$platform[intersect(rownames(study$platform),rownames(targeted_genes)),]
  
  TSS = c()
  for (i in 1:dim(features)[1]){
    if (features$strand[i] == "-"){TSS[i] <- features$tx_end[i]}
    else{TSS[i] <- features$tx_start[i]}
  }
  
  print("Indexing probe by features")
  # params
  ## index meth probes by chr
  pf_chr_colname = "seqnames"
  pf_pos_colname = "start"
  chrs = unique(features[,1])
  chrs_indexed_epic = lapply(chrs, function(chr) {
    print(chr)
    idx = rownames(epic)[epic[[pf_chr_colname]] %in% chr]  
    ret = epic[idx,]
    return(ret)
  })
  names(chrs_indexed_epic) = chrs
  
  ## index probes by gene name
  feat_indexed_probes = epimedtools::monitored_apply(features, 1, function(gene) {
    # gene = features[1,]
    # print(gene)
    chr = gene[[1]]
    meth_platform = chrs_indexed_epic[[chr]]
    ret = dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str) 
    return(ret)
  })
  
  
  features$nb_epic_probes = sapply(feat_indexed_probes[rownames(features)], length)
  features = cbind(features,TSS)
  sub_features =  features[features$nb_epic_probes >= nb_probe_min,]
  
  warning(paste(nrow(features)-nrow(sub_features),"genes were removed due to the followings selection parameters :","\n",
                "window downstream the TSS :",dwn_str,"\n",
                "window upstream the TSS : ", up_str, "\n",
                "Min number of probes :", nb_probe_min, sep = " "))
  
  return(sub_features)
  
}



############################ Visualisation tool for features created by get_features() function##########################################


catalog =  function(selected_features = sub_features){
  
  available_genes = list(catalog="catalog", n_genes = "n_genes")
  
  available_genes[["catalog"]] = noquote(rownames(selected_features))
  available_genes [["n_genes"]] = noquote(paste("Genes selected :", length(available_genes[[1]]),sep = " "))
  warning("Note that quotes were removed for clarity, please select your genes using it for the plotting function.")
  return(available_genes)
  
}

##################### Création des bins #############

plot_selected_probes = function(selected_gene = gene, expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = feat_indexed_probes){

  expr_data = expr_data[,intersect(colnames(expr_data),colnames(meth_data))] #### On ne garde en données que les individus pour lesquels on a méthylation ET 
  
  selected_expr_values = sort(expr_data[selected_gene,]) ##### On sélectionne les données d'expression du gène sélectionné
  selected_probes = intersect(feat_indexed_probes[[selected_gene]],rownames(meth_data)) #### On sélectionne les probes se trouvants sur le gène sélectionné
  tmp = meth_data$platform[selected_probes,]
  
  plot(tmp$Start,rep(1,length(tmp$Start)),yaxt = "n",ylab ="",xlab = "Probes physical coordinates along the gene (in bp) ")
  points(features[selected_gene,"TSS"],1,col="red")
  

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



