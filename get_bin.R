
get_bin_values = function(gene, binlist = list("bin1","bin2","bin3","bin4","bin5","bin6")){

  bin_vals = lapply(binlist, function(bin) {
  tmp_probes = feat_indexed_probes_bin[gene][[bin]]
  tmp_data = meth_lusc$data[intersect(tmp_probes,rownames(meth_lusc$data)),]
  
  return(tmp_data)
  
  })
  
  names(bin_vals) = binlist
  return(bin_vals)
}
  

