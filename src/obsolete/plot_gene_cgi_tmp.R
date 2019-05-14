
plot_gene_cgi <- function(selected_gene = "OTX1",
                          meth_study = meth_lusc,
                          expr_data = trscr_lusc$data,
                          window = c(100000,100000),
                          cgi_platform = cgi_pf,
                          cgi_map = cgi_indexed_probes,
                          features_list = features
                          ){
  
  meth_platform <- meth_study$platform
  meth_data <- meth_study$data
  
  TSS <- features_list[selected_gene,"TSS"]
  
  
  ## get closest cgi + cgi in range
  
  cgi_chr <- cgi_platform[which(as.character(cgi_platform[,"chr"]) == features_list[selected_gene,1]),]
  
  tmp_cgi_of_interest <- cgi_chr[which(abs(cgi_chr[,"center"]-TSS) < mean(window)+5000),]
  cgi_of_interest<- cbind(tmp_cgi_of_interest,abs(tmp_cgi_of_interest[,"center"]-TSS))
  colnames(cgi_of_interest)<-c(colnames(tmp_cgi_of_interest),"dist to TSS")
  
  
  closest_cgi <-cgi_of_interest[which(abs(cgi_of_interest[,"center"]-TSS)==min(abs(cgi_of_interest[,"center"]-TSS),na.rm=T)),]
  
  ## get probes, ordered by physical coordinates
  
  tmp_probes <- unlist(cgi_map[[rownames(closest_cgi)]],use.names = FALSE)
  tmp_sub_epic <- meth_platform[tmp_probes, c(2,3,9)]
  sub_epic<- tmp_sub_epic[order(tmp_sub_epic[,"Start"]),]
  probes <- rownames(sub_epic)
  
  
  ## cols of probes dots
  
  cols<-sapply(rownames(sub_epic), function(probe){
    
    if(sub_epic[probe,"Feature_Type"] == "Island")
    {col = "cornflowerblue"}
    
    
    if(sub_epic[probe,"Feature_Type"] == "S_Shore")
    {col = "chartreuse4"}
    
    
    if(sub_epic[probe,"Feature_Type"] == "S_Shelf")
    {col="gold"}
    
    
    if(sub_epic[probe,"Feature_Type"] == "N_Shore")
    {col = "chocolate"}
    
    
    if(sub_epic[probe,"Feature_Type"] == "N_Shelf")
    {col = "deeppink"}
    
    if(sub_epic[probe,"Feature_Type"] == "opensea")
    {col= "black"}
    
    
    return(col)
  })
  
  layout(matrix(c(4,2,
                  1,3), 2, 2, byrow = TRUE))
  
  ## sorted expr values (pannel 1)
  
  expr_values <- sort(expr_data[selected_gene,intersect(colnames(expr_data),colnames(meth_study$data))])
  
  cols1 <- epimedtools::monitored_apply(t(t(names(expr_values))),1, function(i){
    if(meth_study$exp_grp[i,14] == "tumoral"){col<- "red"}
    else {col<- "blue"}})
  
  tmp = cbind(expr_values,1:length(expr_values))
  expression_plot = plot(tmp,
                         ylab = paste("Patient index ( ordered by level of expression)"),
                         xlab = paste("Expression level"),
                         main = "Expression/Transcription",
                         col = cols1,
                         xlim=c(0,20))
  
  
  
  ## probes location ( pannel 2)
  
  plot(sub_epic$Start, rep(1,length(sub_epic$Start)),
       pch=19,
       xlim=c(min(sub_epic$Start)-50,
              max(sub_epic$End)+50),
       cex=1.1,
       yaxt="n",
       main = paste0("nprobes = ",nrow(sub_epic)),
       xlab="Coordinates (in bp)",
       ylab="",
       ylim=c(0,2),
       col = cols)
  
  text(TSS,0.4,labels = paste0("TSS : ", TSS),cex = 0.7)
  points(TSS, 0.5, pch=9, col="red")
  abline(h=0.5,lty=3)
  
  

  
  if(features[selected_gene,"strand"]== "-"){
    arrows(max(sub_epic$End)-50,1.5,min(sub_epic$Start)+50,1.5)
  }
  
  if(features[selected_gene,"strand"] == "+"){
    arrows(min(sub_epic$Start)+50,1.5,max(sub_epic$End)-50,1.5)
  }
  

  
  
  ## get values for heatmaps ( pannel 3)
  
  
  methvals_of_interest <- meth_data[probes,intersect(colnames(meth_data),names(expr_values))]
  sorted_vals = methvals_of_interest[,order(match(colnames(methvals_of_interest),names(expr_values)))]
  
  colors=c("cyan", "black", "red")
  cols = colorRampPalette(colors)(100)
  
  breaks <- c(seq(0,0.33,length=35),
              seq(0.34,0.66,length=33),
              seq(0.67,1,length=33))
  
  image(sorted_vals,
        axes=FALSE,
        col=cols,
        breaks = breaks,
        ylab="Patients",
        xlab="probes",
        main= "Methylation")
  
  mtext(paste(noquote(selected_gene),"Methylation profile",sep= " "), outer=TRUE,  cex=1.5, line=-2)
  
  ## Legend pannel (4)
  
  plot(1,0,xaxt="n",yaxt="n",xlab="",ylab="",col="white")
  
  legend("right",
         c(paste0("Island : ",sum(sub_epic[,"Feature_Type"]== "Island")),
           paste0("S_Shore : ",sum(sub_epic[,"Feature_Type"]== "S_Shore")),
           paste0("S_Shelf : ",sum(sub_epic[,"Feature_Type"]== "S_Shelf")),
           paste0("N_Shore : ",sum(sub_epic[,"Feature_Type"]== "N_Shore")),
           paste0("N_Shelf : ",sum(sub_epic[,"Feature_Type"]== "N_Shelf")),
           paste0("Opensea : ",sum(sub_epic[,"Feature_Type"]== "opensea"))),
         xpd = TRUE,
         pch=20,
         col=c("cornflowerblue","chartreuse4","gold","chocolate","deeppink","black"))
  
  legend("bottomleft",
         c("tumoral","normal"),
         xpd = TRUE,
         pch=20,
         col=c("red","blue"))
  
  
  
  return(sub_epic)
  
}  
