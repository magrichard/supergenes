#---------------------------------------------
  #'Generate simulated dataset.
  #'
  #'A bedfile (gene information).
  #'
  #' @param n_genes A number of genes to simulate
  #'
  #'@export
  generate_fakestudy <- function(n_genes) {
  start = unique(round(stats::runif(n_genes)*1000000))
  gene_list = data.frame(chr='chr1',
                         start =start,
                         stop = start + round(abs(stats:: rnorm(length(start)) * 10000)),
                         gene_id='gene',
                         score = 'NA',
                         strand = ifelse(stats::runif(length(start))>0.5, '+', '-'),
                         ref_genome='simu',
                         stringsAsFactors=FALSE)
  gene_list$gene_id = paste0('gene_', 1:length(start))
  rownames(gene_list) = gene_list$gene_id
  utils::head(gene_list)
  dim(gene_list)
  return(gene_list)
  }
#---------------------------------------------
  #'Generate simulated dataset.
  #'
  #'A bedfile (gene information).
  #'
  #' @param n_genes A number of genes to simulate
  #'
  #'@export
  generate_fakestudy <- function(n_genes) {
  start = unique(round(stats::runif(n_genes)*1000000))
  gene_list = data.frame(chr='chr1',
                         start =start,
                         stop = start + round(abs(stats:: rnorm(length(start)) * 10000)),
                         gene_id='gene',
                         score = 'NA',
                         strand = ifelse(stats::runif(length(start))>0.5, '+', '-'),
                         ref_genome='simu',
                         stringsAsFactors=FALSE)
  gene_list$gene_id = paste0('gene_', 1:length(start))
  rownames(gene_list) = gene_list$gene_id
  utils::head(gene_list)
  dim(gene_list)
  return(gene_list)
  }
