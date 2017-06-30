``` r
# Make RColorBrewer colors transparent
rcbalpha <- function(alpha=1, ...){
  aa <- data.frame(t(col2rgb(RColorBrewer::brewer.pal(...), alpha=F)))/255
  #aa$alpha <- alpha
  rgb(aa, alpha=alpha)
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


# tetrahedral plot
source('R/dichtcp.R')

# simulate two groups of data
simdich <- function(N=50, sgsqsrate=10, multiplier=c(0.95, 1.05), effsize=NULL){

  musA <- runif(4, 1e-6, 1e0) # vector of means for group A
  musA <- runif(4, 1, 10)
  sgsqs <- rexp(4, sgsqsrate) # vector of standard deviations
  
  groupA <- matrix(NA, nrow=N, ncol=4)
  groupA[,1] <- rlnorm(N, meanlog=log(musA[1]), sdlog=sgsqs[1])
  groupA[,2] <- rlnorm(N, meanlog=log(musA[2]), sdlog=sgsqs[2])
  groupA[,3] <- rlnorm(N, meanlog=log(musA[3]), sdlog=sgsqs[3])
  groupA[,4] <- rlnorm(N, meanlog=log(musA[4]), sdlog=sgsqs[4])
  
  if(is.null(effsize)){
    #musB <- musA*runif(4, .8, 1.2) # vector of means for group B
    musB <- musA * runif(4, multiplier[1], multiplier[2])
  }else{
    # effsize: difference between means in units of standard deviations
    # (mahalanobis distance, approximately)
    sds <- apply(groupA, 2, sd)
    musB <- effsize*sds/2 + musA
  }
  
  groupB <- matrix(NA, nrow=N, ncol=4)
  groupB[,1] <- rlnorm(N, meanlog=log(musB[1]), sdlog=sgsqs[1])
  groupB[,2] <- rlnorm(N, meanlog=log(musB[2]), sdlog=sgsqs[2])
  groupB[,3] <- rlnorm(N, meanlog=log(musB[3]), sdlog=sgsqs[3])
  groupB[,4] <- rlnorm(N, meanlog=log(musB[4]), sdlog=sgsqs[4])
  
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

# testing
#simvals <- runif(100,0,10)
#aa <- lapply(simvals, simdich, N=50, sgsqsrate=10)
#mahd <- unlist(lapply(aa, function(x) sqrt(mahalanobis(colMeans(x[1:50,]), colMeans(x[51:100,]), cov(x[1:50,])) )) )
#plot(simvals~mahd)

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

  grouping <- gsub('[0-9]','', rownames(M))
  
  adonis(M~grouping, ...)
  }

# split data and run volume overlap
voloverlaptest <- function(dat){
  tcsdat <- suppressWarnings(colspace(dat, space='tcs'))
  gA <- tcsdat[1:(dim(dat)[1]/2),]
  gB <- tcsdat[(dim(dat)[1]/2+1):(dim(dat)[1]),]
  
  voloverlap(gA, gB)
}

# split data, get centroids, get color distance
centroidist <- function(dat){
  gA.c <- colMeans(dat[1:(dim(dat)[1]/2),])
  gB.c <- colMeans(dat[(dim(dat)[1]/2+1):(dim(dat)[1]),])
  
  suppressWarnings(coldist(rbind(gA.c, gB.c), achro=FALSE, qcatch='Qi'))$dS
}
```

False positives and power
-------------------------

``` r
effs <- c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3)
timeseach <- 100

effsims <- rep(effs, each=timeseach)

simulatedata <- lapply(effsims,
                       function(x)
                       simdich(N=50, sgsqsrate=10, multiplier=NULL, effsize=x)
                       )


simulatecoldist <- parallel::mclapply(simulatedata, function(x) {
  Y <- suppressWarnings(coldist(x, achro=FALSE, qcatch='Qi'))
  Y$comparison <- NA
  Y$comparison[grepl('A', Y$patch1) & grepl('A', Y$patch2)] <- 'intra.A'
  Y$comparison[grepl('B', Y$patch1) & grepl('B', Y$patch2)] <- 'intra.B'
  Y$comparison[grepl('A', Y$patch1) & grepl('B', Y$patch2)] <- 'inter'
  Y
  }, mc.cores=6)
```

Run stuff

``` r
adonissim <- parallel::mclapply(simulatecoldist, adoniscoldist, mc.cores=6)
vovsim <- parallel::mclapply(simulatedata, voloverlaptest, mc.cores=6)
centdist <- unlist(lapply(simulatedata, centroidist))
gc()
```

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  6643862 354.9    9968622  532.4   9968622  532.4
    ## Vcells 98630284 752.5  160900964 1227.6 131165552 1000.8

Plot stuff

![](../output/figures/final/final_figunnamed-chunk-2-1.png)

color legend:

-   dark colors: methods disagree (BAD)
-   light colors: methods agree (GOOD)

-   light blue: adonis and centroid distance &gt; 1 (GOOD)
-   dark blue: adonis significant, centroid distance &lt; 1 (BAD)
-   dark red: adonis non-significant, centroid distance &gt; 1 (BAD)
-   light red: adonis and centroid distance &lt; 1 (GOOD) ![](../output/figures/final/final_figunnamed-chunk-3-1.png)![](../output/figures/final/final_figunnamed-chunk-3-2.png)![](../output/figures/final/final_figunnamed-chunk-3-3.png)

``` r
sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.4
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] vegan_2.4-3          lattice_0.20-35      permute_0.9-4       
    ## [4] scatterplot3d_0.3-40 pavo_1.1.0          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.11       cluster_2.0.6      knitr_1.15.1      
    ##  [4] magrittr_1.5       maps_3.2.0         magic_1.5-6       
    ##  [7] MASS_7.3-47        geometry_0.3-6     stringr_1.2.0     
    ## [10] tools_3.4.0        parallel_3.4.0     grid_3.4.0        
    ## [13] nlme_3.1-131       mgcv_1.8-17        htmltools_0.3.6   
    ## [16] yaml_2.1.14        rprojroot_1.2      digest_0.6.12     
    ## [19] Matrix_1.2-9       RColorBrewer_1.1-2 mapproj_1.2-5     
    ## [22] codetools_0.2-15   rcdd_1.1-13        evaluate_0.10     
    ## [25] rmarkdown_1.5      stringi_1.1.5      compiler_3.4.0    
    ## [28] backports_1.0.5
