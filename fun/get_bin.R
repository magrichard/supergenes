
########### for loop method ( faster) ----------------------------


get_bin_values1 <- function(data = meth_lusc$data, indexes = feat_indexed_probes_bin, binlist = c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6"), features = features,...) {
  values_per_genes_per_bins <- list()

  for (i in 1:length(binlist)) {
    bin <- binlist[i]

    values_per_genes_per_bins[[paste("bin", i, sep = "")]] <- lapply(rownames(features), function(gene) {
      tmp_probes <- indexes[[gene]][[bin]]
      tmp_data <- data.frame(data[intersect(tmp_probes, rownames(data)), ])

      return(tmp_data)
    })

    names(values_per_genes_per_bins[[i]]) <- rownames(features)
  }

  return(values_per_genes_per_bins)
}



################ sapply method ( slower) ------------------------




get_bin_values2 <- function(data = meth_lusc$data, indexes = feat_indexed_probes_bin, binlist = c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6"), features = features,...) {
  
  values_per_genes_per_bins <- list()
  
  ret = sapply(binlist, function(bin) {
    
      values_per_genes_per_bins[[bin]] <- lapply(rownames(features), function(gene) {
      tmp_probes <- indexes[[gene]][[bin]]
      tmp_data <- data[intersect(tmp_probes, rownames(data)),]
      
      return(tmp_data)
    })
    
    names(values_per_genes_per_bins[[bin]]) <- rownames(features)
    return(values_per_genes_per_bins)
  })
  names(ret) = binlist
  
  return(ret)
}




