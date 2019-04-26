############## Tumoral ################
feat_ind <- get_indexed_binreg()
map_binreg<-reduce_map(feat_ind,c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))


means_per_regions_per_genes_per_patient_t <- reduce_rows(meth_tumoral,map_binreg, mean ,na.rm=T)
means_t <- subset_vals_per_bins(data = meth_tumoral,
                             values_per_patient = means_per_regions_per_genes_per_patient_t,
                             fun = mean,
                             binlist=c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))

meth_heatmap(means_t, main = "mean of means tumoral")


sd_per_regions_per_genes_per_patient_t <- reduce_rows(meth_tumoral,map_binreg, sd ,na.rm=T)
sds_t <- subset_vals_per_bins(data = meth_tumoral,
                              values_per_patient = sd_per_regions_per_genes_per_patient_t,
                              fun = sd,
                              binlist=c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))
meth_heatmap(sds_t, main = "mean of sd tumoral")



rsd_per_regions_per_genes_per_patient_t <- reduce_rows(meth_tumoral,map_binreg, rsd ,na.rm=T)
rsds_t <- subset_vals_per_bins(data = meth_tumoral,
                            values_per_patient = rsd_per_regions_per_genes_per_patient_t,
                            fun = rsd,
                            binlist=c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))
meth_heatmap(rsds_t, main = "mean of rsd tumoral")

############# Healthy

means_per_regions_per_genes_per_patient_h<- reduce_rows(meth_normal,map_binreg, mean ,na.rm=T)
means_h <- subset_vals_per_bins(data = meth_normal,
                              values_per_patient = means_per_regions_per_genes_per_patient_h,
                              fun = mean,
                              binlist=c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))

meth_heatmap(means_h, main = "mean of means healthy")


sd_per_regions_per_genes_per_patient_h <- reduce_rows(meth_normal,map_binreg, sd ,na.rm=T)
sds_h <- subset_vals_per_bins(data = meth_normal,
                            values_per_patient = sd_per_regions_per_genes_per_patient_h,
                            fun = sd,
                            binlist=c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))
meth_heatmap(sds_h, main = "mean of sd healthy")



rsd_per_regions_per_genes_per_patient_h <- reduce_rows(meth_normal,map_binreg, rsd ,na.rm=T)
rsds_h <- subset_vals_per_bins(data = meth_normal,
                             values_per_patient = rsd_per_regions_per_genes_per_patient_h,
                             fun = rsd,
                             binlist=c("bin1","bin2","bin3","INTER","UTR5","INTRON","CDS","UTR3"))
meth_heatmap(rsds_h, main = "mean of rsd healthy")










boxplot_res(means_h,means_t)
boxplot_res(sds_h,sds_t)
boxplot_res(rsds_h,rsds_t)
