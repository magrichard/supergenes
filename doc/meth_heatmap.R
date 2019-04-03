#' Draw a heatmap of results
#' 
#' @description \code{meth_heatmap} draws an heatmap using \code{gplots::heatmap.2} with some specific parameters and with pre-computed color-gradient. 
#'  
#' @param data map-reduced data based, basically the output of \code{subset_vals_per_bins_donc}. 
#' @param ... arguments to be passed in the \code{gplots::heatmap.2} function.
#' @return A heatmap with no dendogram and with color gradient such as green is for low methylation, black for average and red for highly methylated features.
#'@example vars_per_bins_per_genes_per_patient = reduce_rows(meth_lusc$data,binmap,var,na.rm=T)
#' vars = subset_vals_per_bins(values_per_patient = vars_per_bins_per_genes_per_patient)
#' meth_heatmap(vars)
#' @export 



meth_heatmap <- function(data = means, ...){
  
  data = data[order(data[,7]),]
  data = data[,-7]
  dendrogram = "none"
  Rowv = NULL
  Colv = NULL
  
  
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(20)
  
  foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, trace="none", col=cols, main=paste0("Mean of values (", nrow(data), " genes x ", ncol(data), "bins)"), mar=c(10,5), useRaster=TRUE)
  
  return(foo)
}