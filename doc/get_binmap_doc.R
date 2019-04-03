#' Get a map similar to \code{feat_indexed_probes_bins} with only one depth
#' 
#' @description \code{get_binmap} returns a one level depth S3 object which can be used as input in \code{reduce_rows} to map-reduce methylation data based on bins and features 
#'  
#' @param feat_indexed_probes_bin output of \code{get_features_bins}, a S3 object of 2 level depth giving probes per bins for each features
#' @param binlist A Vector object containing names of bins as computed in \code{feat_indexed_probes_bins}
#' @return A one level depth S3 object such as output$gene1.bin1, output$gene1.bin2,...., output$gene2.bin1
#' @example binmap = get_binmap(feat_indexed_probes_bin = feat_indexed_probes_bin, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"))
#' @export 
get_binmap <- function(feat_indexed_probes_bin = feat_indexed_probes_bin, binlist = c("bin1","bin2","bin3","bin4","bin5","bin6"))  {
  bin_indexed_probes = lapply(feat_indexed_probes_bin, function(g){
    g[binlist]
  })
  bin_indexed_probes = unlist(bin_indexed_probes, recursive = FALSE)
}