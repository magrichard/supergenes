###################################  Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################

#' Get probes indexed by features
#' 
#' @description \code{get_features} returns a table in BED format of the features list given and load a \code{feat_indexed_probes} S3 object into global environnment, containing every
#' probes indexed by features.
#' 
#' @param targeted_genes A vector of strings containing names of features of interest.
#' @param study A \code{study} object from which transcription and methylation database are extracted.
#' @param up_str Integer. Length in bp of the upstream TSS-relative egion considered to allocate probes to features.
#' @param dwn_str Integer. Length in bp of the downstream TSS-relative region considered to allocate probes to features.
#' @param nb_probe_min Integer. Minimal number of probes in the considered region to keep feature in the list.
#' @return \code{get_features} return a table of type \code{data.frame} in BED format containing \code{targeted_genes} information.
#' @return \code{feat_indexed_probes} Directly loaded into the global environnment, a S3 object containing EPIC probes indexed by features inputed in \code{target_genes} given
#' \code{up_str}, \code{dwn_str} and \code{nb_probe_min} parameters.
#' @example features = get_features(penda_superdown_deregulated, trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1)
#' @export 
#' 
  get_features <- function(targeted_genes, study, up_str = 2500, dwn_str = 2500, nb_probe_min = 1) { 
    print('Creating a list of features...')
    features <- study$platform[intersect(rownames(study$platform), targeted_genes), ]
    
    TSS <- c()
    for (i in 1:dim(features)[1]) {
      if (features$strand[i] == '-') {
        TSS[i] <- features$tx_end[i]
      }
      else {
        TSS[i] <- features$tx_start[i]
      }
    }
    
    print('Indexing probe by features')

    pf_chr_colname <- 'seqnames'
    pf_pos_colname <- 'start'
    chrs <- unique(features[, 1])
    chrs_indexed_epic <- lapply(chrs, function(chr) {
      print(chr)
      idx <- rownames(epic)[epic[[pf_chr_colname]] %in% chr]
      ret <- epic[idx, ]
      return(ret)
    })
    names(chrs_indexed_epic) <- chrs
    

    feat_indexed_probes <<- epimedtools::monitored_apply(features, 1, function(gene) {

      chr <- gene[[1]]
      meth_platform <- chrs_indexed_epic[[chr]]
      tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
      sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]

      return(tmp_probes)
    })
    
    
    features$nb_epic_probes <- sapply(feat_indexed_probes[rownames(features)], length)
    features <- cbind(features, TSS)
    sub_features <- features[features$nb_epic_probes >= nb_probe_min, ]
    
    return(sub_features) 
  }

  
  

  
  
  

  