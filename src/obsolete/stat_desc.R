rsd_probes <- epimedtools::monitored_apply(meth_lusc$data,1, rsd,na.rm=T)
sd_probes <- epimedtools::monitored_apply(meth_lusc$data,1, sd,na.rm=T)
means_probes <- epimedtools::monitored_apply(meth_lusc$data,1, mean,na.rm=T)
stats_probes <- cbind(sd_probes,means_probes,rsd_probes)
stats_probes<- stats_probes[which(!is.na(stats_probes[,1])),]



genes_values <- reduce_rows(meth_lusc$data,feat_indexed_probes, mean,na.rm=T)
sd_genes <- epimedtools::monitored_apply(genes_values,1, sd,na.rm=T)
means_genes<- epimedtools::monitored_apply(genes_values,1, mean,na.rm=T)
rsd_genes <- epimedtools::monitored_apply(genes_values,1, rsd,na.rm=T)
stats_genes <- cbind(sd_genes,means_genes,rsd_genes)
meth_heatmap(stats_genes)
