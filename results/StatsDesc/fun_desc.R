Plot_Categ_Continu<-function(data,categ,continuous, bins = 50, alpha = .3, binwidth = 0.005){
  require(ggplot2)
  
  hist = ggplot(data, aes(x=continuous, y = ..density.., fill= categ)) + geom_histogram(position = "dodge")+stat_bin(bins=bins,binwidth = binwidth) + xlim(min(continuous),max(continuous))
  dens = ggplot(data, aes(x=continuous, fill=categ)) + geom_density(alpha=alpha)
  box = ggplot(data, aes(x=categ, y=continuous, fill=categ)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  
  res = list("hist" = hist,"density" = dens,"box" = box)
  
  return(res)

  
  
}
