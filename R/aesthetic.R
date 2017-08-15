# Make RColorBrewer colors transparent
rcbalpha <- function(alpha=1, ...){
  aa <- data.frame(t(col2rgb(brewer.pal(...), alpha=F)))/255
  rgb(aa, alpha=alpha)
}

# Define the range for plots
plotrange <- function(x, log=TRUE){
  res <- range(x)
  res[1] <- floor(res[1])
  res[2] <- ceiling(res[2])
  if(log && res[1] == 0)
    res[1] <- range(x)[1]*0.8
  
  if(log && res[2] < 10)
    res[2] <- 10
  
  res
}

# significant/not significant histogram 
yesnohist <- function(x, xlab=""){
  bq <- 0.5  
  while(length(seq(0, ceiling(max(x)),by=bq)) < 10) 
    bq <- bq/2
  
  yes <- hist(x[adonisP], breaks=seq(0, ceiling(max(x)),by=bq), plot=F)
  no <- hist(x[!adonisP], breaks=seq(0, ceiling(max(x)),by=bq), plot=F)
  
  plot(c(0,yes$breaks), c(0,yes$counts,0), type='s',col=palette[1], lwd=5,
       ylab='Frequency', xlab=xlab,
       ylim=range(c(yes$counts,no$counts))* c(0,1.1),
       xlim=c(0, max(c(yes$breaks, no$breaks)*1.1)))
  lines(c(0,no$breaks), c(0,no$counts,0), col=palette[2], type='s', lwd=5)
  
  legend('topright', col=palette[1:2], c('significant  ','not significant  '), 
         pch=15, pt.cex=2.5, bty='n')
}
