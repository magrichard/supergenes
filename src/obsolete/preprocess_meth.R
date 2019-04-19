

process_meth_data<-function(cnv_processed_data = cnv_processed,
                            meth_data = meth_lusc$data,
                            probes_index = feat_indexed_probes,
                            features_list = features){



targeted_cnv <- cnv_processed_data[intersect(rownames(cnv_processed_data),rownames(features_list)),
                                   intersect(colnames(cnv_processed_data),colnames(meth_data))]
  
    methvals_per_genes <- epimedtools::monitored_apply(mod = 10, t(t(names(probes_index))), 1, function(f) {
      
      
    tmp_probes <- intersect(probes_index[[f]],rownames(meth_data))
    selected_meth <- meth_data[tmp_probes,intersect(colnames(targeted_cnv),colnames(meth_data))]  #get data for the f gene
    
    
    
    dummy_matrix <- matrix(nrow = nrow(selected_meth),ncol = ncol(selected_meth),
                       dimnames = list(rownames(selected_meth),colnames(selected_meth)))
    
    
    tmp<-apply(t(t(colnames(selected_meth))),1,function(i) {
      if(targeted_cnv[f,i] != 1){dummy_matrix[,i] <- NA }      #process values per patient for the f gene
      else{dummy_matrix[,i]<-selected_meth[,i]}
    })
    
    
    methvals <- do.call(cbind,tmp)                      #processed meth values for a given f gene
    return(methvals)
    
  })
  
  meth_processed <- do.call(rbind,methvals_per_genes) 
  colnames(meth_processed)<- colnames(targeted_cnv)
   
  return(meth_processed)
}

  

 