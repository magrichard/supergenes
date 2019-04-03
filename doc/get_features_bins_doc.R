
#' Get probes indexed by bins and by features
#' 
#' @description \code{get_features_bins} returns\code{feat_indexed_probes_bins} S3 object into global environnment, containing every
#' probes indexed by features and by bins, with statics bins, even if it might change, such as :
#'  bin1 = -2500bp, -1000bp
#'  bin2 = -1000bp, -500bp
#'  bin3 = -500bp, TSS
#'  bin4 = TSS, +500bp
#'  bin5 = +500bp, +1000bp
#'  bin6 = +1000bp, +2500bp
#'  
#' @param targeted_genes A vector of strings containing names of features of interest.
#' @param study A \code{study} object from which transcription and methylation database are extracted.
#' @param up_str Integer. Length in bp of the upstream TSS-relative egion considered to allocate probes to features.
#' @param dwn_str Integer. Length in bp of the downstream TSS-relative region considered to allocate probes to features.
#' @param nb_probe_min Integer. Minimal number of probes in the considered region to keep feature in the list.
#' @return \code{get_features} return a table of type \code{data.frame} in BED format containing \code{targeted_genes} information.
#' @return \code{feat_indexed_probes} Directly loaded into the global environnment, a S3 object containing EPIC probes indexed by features inputed in \code{target_genes} given
#' \code{up_str}, \code{dwn_str} and \code{nb_probe_min} parameters.
#' @example feat_indexed_probes_bin = get_features_bins(penda_superdown_deregulated,index_pathology = 2,trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1)
#' @export 
get_features_bins <- function(targeted_genes, study = trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1, ...) { #### In my uses-case, 1 for LUAD & 2 for LUSC
  print("Creating a list of features...")
  features <- study$platform[intersect(rownames(study$platform), targeted_genes), ]
  
  
  
  
  print("Indexing probe by features and by bins...")

  pf_chr_colname <- "seqnames"
  pf_pos_colname <- "start"
  chrs <- unique(features[, 1])
  chrs_indexed_epic <- lapply(chrs, function(chr) {
    print(chr)
    idx <- rownames(epic)[epic[[pf_chr_colname]] %in% chr]
    ret <- epic[idx, ]
    return(ret)
  })
  names(chrs_indexed_epic) <- chrs
  
  
  

  feat_indexed_probes_bin <- epimedtools::monitored_apply(features, 1, function(gene) {

    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[chr]]
    tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) - 2500
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) + 1000
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 1500
    bin1 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) - 1000
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) + 500
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin2 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) - 500
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) + 0
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin3 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) + 0
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) - 500
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin4 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) + 500
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) - 1000
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 500
    bin5 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    tmp_gene <- gene
    tmp_gene[[6]] <- "+"
    if (gene[[6]] == "+") {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[2]]) + 1000
    } else {
      tmp_gene[[2]] <- as.numeric(tmp_gene[[3]]) - 2500
    }
    tmp_u_bin <- 0
    tmp_d_bin <- 1500
    bin6 <- dmprocr::get_probe_names(tmp_gene, sub_epic, pf_chr_colname, pf_pos_colname, tmp_u_bin, tmp_d_bin)
    
    
    ret <- list(sub_epic = sub_epic, bin1 = bin1, bin2 = bin2, bin3 = bin3, bin4 = bin4, bin5 = bin5, bin6 = bin6)
    return(ret)
    
    
  })
  
  
  
  
  return(feat_indexed_probes_bin)
}