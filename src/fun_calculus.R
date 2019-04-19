

################## relative standard deviation ##############


rsd <- function(x,...) {
  
  sd(x,...) / abs(mean(x,...))
  
}
