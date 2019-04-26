
features <- get_features(penda_superconserved, study = trscr_lusc, up_str = 7500, dwn_str = 7500)
layout(matrix(1:2,1), respect=TRUE)

plot_mean_vs_sd <- function(genes_list = feat_indexed_probes,
                            data = meth_tumoral,...){
  

  
  poi <-  unique(unlist(genes_list))
  meth <- data[intersect(rownames(data),poi),]
  
  rsd_probes <- epimedtools::monitored_apply(meth,1, rsd,na.rm=T)
  sd_probes <- epimedtools::monitored_apply(meth,1, sd,na.rm=T)
  means_probes <- epimedtools::monitored_apply(meth,1, mean,na.rm=T)
  
  tmp_stats <- cbind(means_probes[!is.na(means_probes)],sd_probes[!is.na(sd_probes)],rsd_probes[!is.na(rsd_probes)])
  stats <- tmp_stats[order(tmp_stats[,1]),]
  colnames(stats)<- c("means","sd","rsd")
  
  
  plot(stats[,"means"],stats[,"sd"],
       xlab = paste0("mean per probe (nprobes = ",nrow(stats)," + ",nrow(meth)-nrow(stats)," NAs )"),
       ylab="sd per probe",...)
  
  return(stats)
}



