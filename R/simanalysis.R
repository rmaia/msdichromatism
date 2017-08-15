
# make distance matrix and run adonis
adoniscoldist <- function(x, ...){
  coldistres <- as.matrix(rbind(x[ ,c(1,2,3)], x[ ,c(2,1,3)]))
  uniquepatches <-  unique(c(coldistres[,1], coldistres[,2]))
  
  M <- matrix(nrow=length(uniquepatches), ncol=length(uniquepatches))
  
  rownames(M) <- colnames(M) <- uniquepatches
  
  M[coldistres[,1:2] ] <- coldistres[,3]
  M[coldistres[,2:1] ] <- coldistres[,3]
  
  class(M) <- 'numeric'
  M[is.na(M)] <- 0
  
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
