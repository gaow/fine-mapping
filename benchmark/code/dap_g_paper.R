## args = commandArgs(trailingOnly=TRUE)
## index = args[1]
## X = as.matrix(read.table("sim.r.geno"))
n = dim(X)[1]      # number of observations
p = dim(X)[2]      # number of variables
freq = 0.05

region = rbinom(91,1,freq)
amplitude = 0.6

get.beta<-function(v){
  rstv = rep(0,11)
  if(v==1){
    index = sample(1:11,1)
    rstv[index]=rnorm(1,sd=amplitude);

  }
  return(rstv)
}

beta = as.vector(sapply(region, get.beta))

y.sample <- function() X %*% beta + rnorm(n)

Y = matrix(ncol=1,y.sample())
B = matrix(ncol=1,beta)
