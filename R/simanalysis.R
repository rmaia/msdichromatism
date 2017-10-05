
# make distance matrix and run adonis
adoniscoldist <- function(x, ...){
  M <- coldist2mat(x)[['dS']]  
  grouping <- as.factor(gsub('[0-9]','', rownames(M)))
  
  M <- as.dist(M)
  
  adonis(M~grouping, ...)
}

# split data and run volume overlap
voloverlaptest <- function(dat, jnd2xyzres=FALSE){
  if(jnd2xyzres){
    tcsdat <- dat
  }else{
    tcsdat <- suppressWarnings(colspace(dat, space='tcs'))
  }
  gA <- tcsdat[1:(dim(dat)[1]/2),]
  gB <- tcsdat[(dim(dat)[1]/2+1):(dim(dat)[1]),]
  
  voloverlap(gA, gB)
}

# geometric mean
gmean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

# split data, get centroids based on geometric means, get color distance
centroidist <-  function(dat){
  gA.c <- apply(dat[1:(dim(dat)[1]/2),], 2, gmean) 
  gB.c <- apply(dat[(dim(dat)[1]/2+1):(dim(dat)[1]),], 2, gmean) 
  suppressWarnings(coldist(rbind(gA.c, gB.c), achro=FALSE, qcatch='Qi'))$dS
}
