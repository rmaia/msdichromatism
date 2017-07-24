# simulate lognormal with realized mean and standard deviation on the log scale
actuallognorm <- function(N, amean, asd){
  location <- log(amean^2 / sqrt(asd^2 + amean^2))
  shape <- sqrt(log(1 + (asd^2 / amean^2)))
  rlnorm(n=N, meanlog=location, sdlog=shape)
}

# simulate two groups of data
simdich <- function(N=50, sgsqsrate=1, sdmeanratio=TRUE, sdmeanratiovalues=c(0,0.5), multiplier=c(0.95, 1.05), effsize=NULL){
  
  musA <- runif(4, 1, 10) # vector of means for group A
  
  if(sdmeanratio){
    sgsqs <- runif(4,sdmeanratiovalues[1],sdmeanratiovalues[2])*musA
  }else{
    sgsqs <- rexp(4, sgsqsrate) # vector of standard deviations
  }
  
  groupA <- matrix(NA, nrow=N, ncol=4)
  groupA[,1] <- actuallognorm(N, amean=musA[1], asd=sgsqs[1])
  groupA[,2] <- actuallognorm(N, amean=musA[2], asd=sgsqs[2])
  groupA[,3] <- actuallognorm(N, amean=musA[3], asd=sgsqs[3])
  groupA[,4] <- actuallognorm(N, amean=musA[4], asd=sgsqs[4])
  
  if(is.null(effsize)){
    musB <- musA * runif(4, multiplier[1], multiplier[2])
  }else{
    # effsize: difference between means in units of standard deviations
    # (mahalanobis distance, approximately)
    # up to 4 dimensions (K), so divide the multivariate effect size by sqrt(K)
    sds <- apply(groupA, 2, sd)
    musB <- (effsize/sqrt(4)*sds) + musA
  }
  
  groupB <- matrix(NA, nrow=N, ncol=4)
  groupB[,1] <- actuallognorm(N, amean=musB[1], asd=sgsqs[1])
  groupB[,2] <- actuallognorm(N, amean=musB[2], asd=sgsqs[2])
  groupB[,3] <- actuallognorm(N, amean=musB[3], asd=sgsqs[3])
  groupB[,4] <- actuallognorm(N, amean=musB[4], asd=sgsqs[4])
  
  combined <- data.frame(rbind(groupA,groupB))
  
  colnames(combined) <- c('u','s','m', 'l')
  rownames(combined) <- paste(rep(c('gA','gB'), each=N),1:N, sep='')
  
  attr(combined, 'relative') <- FALSE
  
  simpars <- data.frame(rbind(musA, musB, sgsqs))
  colnames(simpars) <- c('u','s','m', 'l')
  rownames(simpars) <- c('muA','muB','ssq')
  attr(combined, 'simpar') <- simpars
  
  combined
}
