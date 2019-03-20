Plot_Categ_Continu<-function(data,categ,continuous, bins = 50, alpha = .3, binwidth = 0.005, main = "Categ-Continuous", main1 = "Hist", main2 = "dens", main3 = "box"){
  require(ggplot2)
  require(gridExtra)
  
  hist = ggplot(data, aes(x=continuous, y = ..density.., fill= categ)) + geom_histogram(position = "dodge")+stat_bin(bins=bins,binwidth = binwidth) + xlim(min(continuous),max(continuous)) + ggtitle(main1)
  dens = ggplot(data, aes(x=continuous, fill=categ)) + geom_density(alpha=alpha) + ggtitle(main2)
  box = ggplot(data, aes(x=categ, y=continuous, fill=categ)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2) + ggtitle(main3)
  
  
  
  grid.arrange(hist, dens, box, nrow=1, ncol=3, top = main)
  
  
  res = list("hist" = hist,"density" = dens,"box" = box)
  return(res)

  
  
}
