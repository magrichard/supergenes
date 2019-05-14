variability_plot<- function(targeted_genes = penda_superup_deregulated,
                            meth_tumoral_data = meth_tumoral,
                            meth_healthy_data = meth_normal,
                            meth_platform = meth_lusc$platform,
                            window = c(2500,2500),
                            ...)
{
  
  
  ## Data manipulation
  
  features <- get_features(targeted_genes, study = trscr_lusc, up_str = window[2], dwn_str = window[1])
  meth_tumoral_no_sex_chr <- meth_tumoral_data[intersect(rownames(meth_platform[which(meth_platform[,1] != "chrX"
                                                                                      & meth_platform[,1] != "chrY"),]),
                                                         rownames(meth_tumoral_data)),]
  meth_healthy_no_sex_chr <- meth_healthy_data[intersect(rownames(meth_platform[which(meth_platform[,1] != "chrX"
                                                                                      & meth_platform[,1] != "chrY"),]),
                                                         rownames(meth_healthy_data)),]
  
  healthy_of_interest <- meth_healthy_no_sex_chr[intersect(unique(unlist(feat_indexed_probes)),rownames(meth_healthy_no_sex_chr)),]
  tumoral_of_interest <- meth_tumoral_no_sex_chr[intersect(unique(unlist(feat_indexed_probes)),rownames(meth_tumoral_no_sex_chr)),]
  
  ## Get values
  
  sd_tumoral <-apply(tumoral_of_interest,1,sd,na.rm=T)
  mean_tumoral <- apply(tumoral_of_interest,1,mean,na.rm=T)
  rsd_tumoral <- apply(tumoral_of_interest,1,rsd,na.rm=T)
  
  
  sd_healthy<- apply(healthy_of_interest,1,sd,na.rm=T)
  mean_healthy<-apply(healthy_of_interest,1,mean,na.rm=T)
  rsd_healthy <- apply(healthy_of_interest,1,rsd,na.rm=T)
  
  ## Plot
  
  plot(mean_healthy,
       sd_healthy,
       col="blue",
       ylim=c(0,max(c(max(sd_healthy,na.rm=T),max(sd_tumoral,na.rm=T)))),
       xlab = paste0("mean per probe","\n","(n probes = ",length(sd_healthy),"  window : ","[",window[1],",",window[2],"]"," )"),
       ylab = "sd per probe",
       ...)
  
  points(mean_tumoral,sd_tumoral,col="red")
  abline(v = c(mean(mean_healthy,na.rm=T),mean(mean_tumoral,na.rm=T)),col=c("blue","red"),lty = 3)
  abline(h = c(mean(sd_healthy,na.rm=T),mean(sd_tumoral,na.rm=T)),col=c("blue","red"),lty = 3)
  
  arrows(mean_healthy,sd_healthy,mean_tumoral,sd_tumoral,
         length=0.1,
         lwd=0.3)
  
  ## get a table to return for further investigation 
  
  ret <- cbind(mean_healthy,sd_healthy,rsd_healthy,mean_tumoral,sd_tumoral,rsd_tumoral)
  
  return(ret)
  
  
}


