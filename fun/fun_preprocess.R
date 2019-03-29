###################################  Getting probes features for targeted type of genes ( superdown, superup...) and type of lung cancer #######################


get_features <- function(targeted_genes, index_pathology = 2, study, up_str = 2500, dwn_str = 2500, nb_probe_min = 1, ...) { #### In my uses-case, 1 for LUAD & 2 for LUSC
  print("Creating a list of features...")
  targeted_genes <- targeted_genes[targeted_genes[, index_pathology] == 1, ]
  features <- study$platform[intersect(rownames(study$platform), rownames(targeted_genes)), ]

  TSS <- c()
  for (i in 1:dim(features)[1]) {
    if (features$strand[i] == "-") {
      TSS[i] <- features$tx_end[i]
    }
    else {
      TSS[i] <- features$tx_start[i]
    }
  }

  print("Indexing probe by features")
  # params
  ## index meth probes by chr
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

  ## index probes by gene name
  feat_indexed_probes <- epimedtools::monitored_apply(features, 1, function(gene) {
    # gene = features[1,]
    # print(gene)
    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[chr]]
    tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]
    # here compute bin 1 to 6 using sub_epic
    # do_biens
    # ret = list(sub_epic=sub_epic, bin1=..., bin2=..., bin3=..., bin4=..., bin5=..., bin6=...)
    return(tmp_probes)
  })


  features$nb_epic_probes <- sapply(feat_indexed_probes[rownames(features)], length)
  features <- cbind(features, TSS)
  sub_features <- features[features$nb_epic_probes >= nb_probe_min, ]

  warning(paste("\n", nrow(features) - nrow(sub_features), "genes were removed due to the followings selection parameters :", "\n",
    "window downstream the TSS :", dwn_str, "\n",
    "window upstream the TSS : ", up_str, "\n",
    "Min number of probes :", nb_probe_min,
    sep = " "
  ))

  return(sub_features)
}


############################### get_features_bins : same as prior + indexing probes per bins #####################

get_features_bins <- function(targeted_genes, index_pathology = 2, study = trscr_lusc, up_str = 2500, dwn_str = 2500, nb_probe_min = 1, ...) { #### In my uses-case, 1 for LUAD & 2 for LUSC
  print("Creating a list of features...")
  targeted_genes <- targeted_genes[targeted_genes[, index_pathology] == 1, ]
  features <- study$platform[intersect(rownames(study$platform), rownames(targeted_genes)), ]




  print("Indexing probe by features")
  # params
  ## index meth probes by chr
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



  ## index probes by gene name
  feat_indexed_probes_bin <<- epimedtools::monitored_apply(features, 1, function(gene) {
    # gene = features[1,]
    # print(gene)
    chr <- gene[[1]]
    meth_platform <- chrs_indexed_epic[[chr]]
    tmp_probes <- dmprocr::get_probe_names(gene, meth_platform, pf_chr_colname, pf_pos_colname, up_str, dwn_str)
    sub_epic <- meth_platform[tmp_probes, c(pf_chr_colname, pf_pos_colname)]
    # here compute bin 1 to 6 using sub_epic
    # do_bins()
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


  features$nb_epic_probes <- sapply(feat_indexed_probes_bin[rownames(features)], length)
  sub_features <- features[features$nb_epic_probes >= nb_probe_min, ]

  print(paste(nrow(features) - nrow(sub_features), "genes were removed due to the followings selection parameters :",
    "window downstream the TSS :", dwn_str,
    "window upstream the TSS : ", up_str,
    "Min number of probes :", nb_probe_min,
    sep = " "
  ))


  return(sub_features)
}


############## reduce_rows #################

reduce_rows <- function(tmp_meth_data, map, indicator_func2 = mean, ...) {
  dim(tmp_meth_data)
  meth_by_tissues_by_feature <- epimedtools::monitored_apply(mod = 10, t(t(names(map))), 1, function(f) {
    # print(f)
    # f = features[1,]
    # f = features["HACD4",]
    probe_idx <- intersect(rownames(tmp_meth_data), map[[f]])
    if (length(probe_idx) == 0) {
      # print(f)
      tmp_meth_by_tissues_by_feature <- rep(NA, ncol(tmp_meth_data))
      names(tmp_meth_by_tissues_by_feature) <- colnames(tmp_meth_data)
    } else if (length(probe_idx) > 1) {
      tmp_meth_by_tissues_by_feature <- apply(tmp_meth_data[probe_idx, ], 2, indicator_func2, ...)
    } else {
      tmp_meth_by_tissues_by_feature <- sapply(tmp_meth_data[probe_idx, ], indicator_func2, ...)
    }
    return(tmp_meth_by_tissues_by_feature)
  })
  colnames(meth_by_tissues_by_feature) <- names(map)
  meth_by_tissues_by_feature <- t(meth_by_tissues_by_feature)
  return(meth_by_tissues_by_feature)
}








# Genes selected catalog function ( no calculus only visualisation)


catalog <- function(selected_features = features) {
  available_genes <- list(catalog = "catalog", n_genes = "n_genes")

  available_genes[["catalog"]] <- noquote(rownames(selected_features))
  available_genes [["n_genes"]] <- noquote(paste("Genes selected :", length(available_genes[[1]]), sep = " "))
  warning("Note that quotes were removed for clarity, please select your genes using it for the plotting function.")
  return(available_genes)
}



###########################  Creation of bins ######################

do_bins <- function(selected_gene = gene, nbins = 5, ...) {
  window <- c(features$features[selected_gene, "TSS"] - dwn_str, features$features[selected_gene, "TSS"] + up_str)
  window_width <- window[2] - window[1]
  bw <- window_width / nbins


  tmp <- meth_lusc$platform[selected_probes, ]
  bins <- cut(tmp$Start, nbins)
  levels(bins) <- seq(1, nbins, 1)

  bins_coordinates <- matrix(ncol = 2, nrow = nbins)
  for (i in 2:(window_width / bw)) {
    bins_coordinates[1, ] <- c(window[1], window[1] + bw)
    bins_coordinates[i, 1] <- bins_coordinates[i - 1, 2]
    bins_coordinates[i, 2] <- bins_coordinates[i, 1] + bw
  }
  colnames(bins_coordinates) <- c("Start", "End")
  rownames(bins_coordinates) <- levels(bins)

  desc <- cbind(bins, rowMeans(methvals_of_interest, na.rm = TRUE), apply(methvals_of_interest, 1, var, na.rm = TRUE), apply(methvals_of_interest, 1, sd, na.rm = TRUE))
  desc <- rbind(tapply(desc[, 2], desc[, 1], mean), tapply(desc[, 3], desc[, 1], mean), tapply(desc[, 4], desc[, 1], mean))
  rownames(desc) <- c("mean", "Var", "sd")

  res <- list("bins" = bins, "table" = table(bins), "bins_coordinates" = bins_coordinates, "stats" = desc)
  return(res)
}

create_map <- function(indexes = feat_indexed_probes_bin, binlist = c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6"), genes = features) {
  map <- c()

  ret <- sapply(binlist, function(bin) {
      lapply(rownames(features), function(gene) {
      map[[paste(gene, "_", bin, sep = "")]] <- feat_indexed_probes_bin[[gene]][[bin]]
      tmp_names <<- paste(gene, "_", bin, sep = "")
      return(map)
    })

  })
  return(ret)
}
