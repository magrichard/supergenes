#' Subset map-reduced data based on a function such as mean.
#' 
#' @description \code{subset_vals_per_bins} returns a table with bins as columns and features as rows with an additionnal column called "overall" 
#' giving the aggregated metric over the whole dataset
#'  
#' @param values_per_patient map-reduced data based, basically the output of \code{reduce_rows}. 
#' @param binlist A Vector object containing names of bins as computed in \code{feat_indexed_probes_bins}.
#' @param fun the function to apply to aggregate values per patient.
#' @param ... arguments the user wants to pass into the \code{fun} object.
#' @return A table of type data.frame of dimension features*binlist+1 with colnames like c(binlist,"overall")
#' @example means = subset_vals_per_bins(values_per_patient = means_per_bins_per_genes_per_patient, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"), fun = mean, na.rm=T)
#' @export 
subset_vals_per_bins <- function(values_per_patient = means_per_bins_per_genes_per_patient, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"), fun = mean, ...){
  values_per_bins <- apply(values_per_patient,1,fun, na.rm=T) 
  ret <- list()
  
  vals <- sapply(binlist, function(bin){
    ret[[bin]] <- values_per_bins[endsWith(names(values_per_bins),bin)]    
    return(ret)
    
  }
  )
  
  
  tmp = reduce_rows(meth_lusc$data, feat_indexed_probes, mean, na.rm=T) ### note that you need feat_indexed_probes loaded into your environnment
  overall = apply(tmp,1,fun,na.rm=T)
  
  
  vals_per_genes <- do.call(cbind,vals)
  vals_per_genes <- cbind(vals_per_genes,overall)
  colnames(vals_per_genes) <-c(binlist,"overall")
  rownames(vals_per_genes) <- names(feat_indexed_probes_bin) 
  
  return(vals_per_genes)
  
}
