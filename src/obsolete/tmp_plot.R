plot_gene_dmr<- function(selected_gene = "OTX1",
                         expr_data = trscr_lusc$data,
                         meth_study = meth_lusc,
                         dmrs_table = DMR_table,
                         DMR_map  = DMRs100k){
  
  
  
  tmp_DMR <- dmrs_table[names(DMR_map[[selected_gene]]),c("chr","start","end","is.hyper")]
  DMRs_of_interest<- tmp_DMR[order(tmp_DMR[,"start"])]
  
  if(nrow(DMRs_of_interest) == 0) {stop(paste0("There is no DMRs indexed for ",selected_gene,". Check your DMR_map object"))}
  
  probes_on_chr <- epic[which(as.character(epic[,1]) == DMRs_of_interest[["chr"]][1]),]
  
  
  expr_values <- sort(expr_data[selected_gene,intersect(colnames(expr_data),colnames(meth_study$data))])
  meth_values <- meth_study$data[,intersect(colnames(expr_data),colnames(meth_study$data))]
  
  
  
  
  
  
  
  
  
  tmp_DMR <- dmrs_table[names(DMR_map[[selected_gene]]),c("chr","start","end","is.hyper")]
  DMRs_of_interest<- tmp_DMR[order(tmp_DMR[,"start"]),]
  probes_on_chr <- epic[which(as.character(epic[,1]) == DMRs_of_interest[["chr"]][1]),]
  
  print("Indexing methylation values per DMR...")
  vals_per_dmr<-epimedtools::monitored_apply(t(t(rownames(DMRs_of_interest))),1, function(i){
    
    DMR<-DMRs_of_interest[i,]   
    
    
    selected_probes <- probes_on_chr[which(probes_on_chr[,"start"] >= DMR[["start"]] & probes_on_chr[,"end"] <= DMR[["end"]]),]
    selected_meth_values <- meth_values[intersect(rownames(selected_probes),rownames(meth_values)),]
    
    
    if(class(selected_meth_values) != "numeric"){
      DMR_vals <- apply(selected_meth_values,2,fun,na.rm=T)}         #### Check if there is more than 2 selected probes
    
    else{DMR_vals <- rep(NA,length(selected_meth_values))}
    
    
    
    
    #### Else, returning raw values if the function is mean, else return a vector of NAs
    
    
    return(DMR_vals)
    
    
  })
  
  rownames(vals_per_dmr)<- colnames(meth_values)
  sorted_vals = vals_per_dmr[order(match(rownames(vals_per_dmr),names(selected_expr_values))),]
  colnames(sorted_vals)<- rownames(DMRs_of_interest)
  
  cols1 <- epimedtools::monitored_apply(t(t(names(expr_values))),1, function(i){
    if(meth_study$exp_grp[i,14] == "tumoral"){col[i]<- "red"}
    else {col[i]<- "blue"}})
  
  print("Plotting...")
  
  layout(matrix(1:2,1), respect=TRUE)
  
  tmp = cbind(selected_expr_values,1:length(selected_expr_values))
  expression_plot = plot(tmp,
                         ylab = paste("Patient index ( ordered by level of expression)"),
                         xlab = paste("Expression level"),
                         main = "Expression/Transcription",
                         col = cols1,
                         xlim=c(0,20))
  legend("topleft",
         c("tumoral","normal"),
         xpd = TRUE,
         pch=20,
         col=c("red","blue"))
  
  
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(100)
  breaks <- c(seq(0,0.33,length=35),seq(0.34,0.66,length=33),seq(0.67,1,length=33))
  
  image(t(sorted_vals),
        col=cols,
        breaks = breaks,
        axes=FALSE,
        main = paste0("# of DMRs : ",ncol(vals_per_dmr)))
  
  
  return(sorted_vals)
  
  
  
}



  
  