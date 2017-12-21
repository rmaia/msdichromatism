Simulations - power
================

-   [Type I error and and power](#type-i-error-and-and-power)
    -   [Running Analysis](#running-analysis)
    -   [Visualizing Results](#visualizing-results)

``` r
require(pavo)
```

    ## Loading required package: pavo

``` r
require(vegan)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.4-3

``` r
require(RColorBrewer)
```

    ## Loading required package: RColorBrewer

``` r
# load aesthetic functions (plot, make colors transparent)
source('R/aesthetic.R')

# load simulation and analysis functions
source('R/simfoos.R')
source('R/simanalysis.R')
source('R/pausemcl.R')
source('R/mahalanobis.R')
```

Type I error and and power
==========================

We will simulate data with varying effect sizes (centroid distance relative to the variance-covariance matrix, or Mahalanobis distance):

``` r
effs <- c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3)
timeseach <- 200
simN <- 50

effsims <- rep(effs, each=timeseach)

simulatedata <- lapply(effsims,
                       function(x)
                       simdich(N=simN, sgsqsrate=0.5, multiplier=NULL, effsize=x)
                       )

simulatecoldist <- pausemcl(simulatedata, function(x) {
  Y <- suppressWarnings(coldist(x, achro=FALSE, qcatch='Qi'))
  Y$comparison <- NA
  Y$comparison[grepl('A', Y$patch1) & grepl('A', Y$patch2)] <- 'intra.A'
  Y$comparison[grepl('B', Y$patch1) & grepl('B', Y$patch2)] <- 'intra.B'
  Y$comparison[grepl('A', Y$patch1) & grepl('B', Y$patch2)] <- 'inter'
  Y
  } )
```

Validating simulations:

![](../output/figures/final_power_fig_unnamed-chunk-1-1.jpeg)

Verifying that values obtained in the simulation (empirical) are close to what we wanted to simulate (simulated) for the four cones (violet, blue, green, red)

Running Analysis
----------------

Run adonis, volume overlap, calculate distance between centroids:

``` r
gc(verbose=FALSE)
adonissim <- pausemcl(simulatecoldist, adoniscoldist)
vovsim <- pausemcl(simulatedata, voloverlaptest)
centdist <- unlist(lapply(simulatedata, centroidist))
gc(verbose=FALSE)
```

Run Pyke MANOVA:

``` r
# remove grouping variable, reattach attributes
scd2 <- lapply(simulatecoldist,'[', ,1:3, drop=FALSE)
for(i in 1:length(scd2)){
  attributes(scd2[[i]])[
    grep('name', names(attributes(simulatecoldist[[i]])), invert=TRUE, value=TRUE)] <-
    attributes(simulatecoldist[[i]])[
    grep('name', names(attributes(simulatecoldist[[i]])), invert=TRUE, value=TRUE)]
  }

pykesim <- parallel::mclapply(scd2, jnd2xyz, rotate=FALSE, mc.cores=4)
pykelm <-  parallel::mclapply(pykesim, function(x) 
  lm(as.matrix(x) ~ rep(c('gA','gB'), each=50)), mc.cores=4)
pykemanova <- parallel::mclapply(pykelm, function(x) summary(manova(x)), mc.cores=4)

vovpyke <- pausemcl(pykesim, function(x)
  voloverlap(x[1:simN,], x[(simN+1):(simN*2), ]) )
```

Calculate statistics of interest:

``` r
# Pyke MANOVA: P-values, which significant
manovaPval <- unlist(lapply(pykemanova, function(x) x$stats[1,'Pr(>F)']))
manovaP <- manovaPval < 0.05

# Adonis: P-values, which significant, R-squared
adonisPval <- unlist(lapply(adonissim, function(x) x$aov.tab$'Pr(>F)'[1]))
adonisP <- adonisPval < 0.05
adonisR2 <- unlist(lapply(adonissim, function(x) x$aov.tab$'R2'[1])) * 100

# Disagreement in significance between Pyke and Adonis
disagreement <- !manovaP == adonisP

# Proportion of tests significant by effect size, according to each method
pvsimsadonis <- tapply(adonisP, effsims, mean)
pvsimsmanova <- tapply(manovaP, effsims, mean)

# Volume overlap
overlap <- unlist(lapply(vovsim, '[','vboth')) * 100
overlapyke <- unlist(lapply(vovpyke, '[','vboth')) * 100


# Centroids: intra-group distances, inter-group distances, which centroid distances > 1
gmeans <- lapply(simulatecoldist, function(x) tapply(x$dS, x$comparison, mean))
gmeans <- do.call(rbind, gmeans)

intradist <- rowMeans(gmeans[, c('intra.A', 'intra.B')])
interdist <- gmeans[,"inter"]

centroidP <- centdist > 1

# Mahalanobis distance between centroids
mahd <- unlist(lapply(simulatedata, mahalanobis, groups=rep(c("A","B"), each=50), mve=FALSE))

# Color palette for plots
palette <- rcbalpha(0.6, 4, 'Set1')

#sigpal <- as.character(factor(adonisP, labels=palette[1:2]))
sigpal <- as.character(factor(paste(adonisP, centroidP), 
          levels=c("FALSE FALSE", "FALSE TRUE", "TRUE FALSE", "TRUE TRUE"), 
          labels=rcbalpha(0.8, 6, 'RdBu')[c(3,1,6,4)] ))
```

Visualizing Results
-------------------

![](../output/figures/final_power_fig_unnamed-chunk-3-1.jpeg)

The simulation was successful in producing samples that had the desired mahalanobis distance. There is some spread because of the small sample size relative to the dimensionality of the dataset, and for that reason the distances between groups asymptotes before zero.

![](../output/figures/final_power_fig_unnamed-chunk-4-1.jpeg)

Both tests had similar power. For very low effect sizes, Type-I Error rate is close to the desired 0.05 (dashed line).

Both tests are quite sensitive too, with signifcant results when the distance between centroids is of about the same magnitude as the pooled standard deviations. Further, Pyke-MANOVA seems less conservative overall than Adonis, though they are very close to each other at very small effect sizes.

![](../output/figures/final_power_fig_unnamed-chunk-5-1.jpeg)

However, there is some discrepancy in test results. There doesn't seem to be a bias - results are centered around the 1:1 line, difference in P-values from the tests is centered and mostly symmetric around 0. But there are occasions in which results are significant for one test but not the other (red; the space between the dashed lines in the first plot).

So tests have similar power but disagree as to the outcome in terms of what is significant:

    ##        manovaP
    ## adonisP  FALSE   TRUE
    ##   FALSE 0.4910 0.0825
    ##   TRUE  0.0190 0.4075

About 10% divergence in results, maybe not worth worrying about. Note that most of the discrepancy comes from results that are significant in MANOVA but not in Adonis, suggesting again that MANOVA approach is less conservative.

There is disagreement particularly when the effect size is marginal: ![](../output/figures/final_power_fig_unnamed-chunk-7-1.jpeg)

![](../output/figures/final_power_fig_unnamed-chunk-8-1.jpeg)

Even though centroid distance increases with effect size, there's a lot of spread, which is the core of the problem we're trying to address - you can have a huge centroid distance with a small Mahalanobis distance (small separation bewtwen the groups). Note centroid is in log scale.

![](../output/figures/final_power_fig_unnamed-chunk-9-1.jpeg)

R2 increases with increasing effect size, which is good. We can also see that even though a lot of the simulations have a distance between centroids greater than 1 (0.48), they are still not significant (red) according to either approach. Transition from non-siginificant to significant occurs for Mahalanobis Distance between 0.5 and 1.

![](../output/figures/final_power_fig_unnamed-chunk-10-1.jpeg)

This just shows the Pyke transformation is working and that the Euclidean distance between the centroids calculated in this transformed space is identical to the distance between the centroids in JNDs.

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.2
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
    ## [1] RColorBrewer_1.1-2 vegan_2.4-3        lattice_0.20-35   
    ## [4] permute_0.9-4      pavo_1.3.1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.13     cluster_2.0.6    knitr_1.16       magrittr_1.5    
    ##  [5] MASS_7.3-47      maps_3.2.0       magic_1.5-6      geometry_0.3-6  
    ##  [9] stringr_1.2.0    globals_0.10.2   tools_3.4.2      grid_3.4.2      
    ## [13] parallel_3.4.2   nlme_3.1-131     mgcv_1.8-20      htmltools_0.3.6 
    ## [17] yaml_2.1.14      rprojroot_1.2    digest_0.6.12    Matrix_1.2-11   
    ## [21] pbmcapply_1.2.4  mapproj_1.2-5    codetools_0.2-15 rcdd_1.2        
    ## [25] evaluate_0.10.1  rmarkdown_1.6    stringi_1.1.5    compiler_3.4.2  
    ## [29] backports_1.1.0  future_1.6.1     listenv_0.6.0

plots for publication:

``` r
############
# EXAMPLES #
############
pdf(height=4*1.3, width=7*1.3, file='figures/exampletetra.pdf')
par(mfrow=c(1,2))
eg1 <- simulatedata[[10]]
attr(eg1, 'conenumb') <- '4'
attr(eg1, 'relative') <- FALSE
attr(eg1, 'qcatch') <- 'Qi'
class(eg1) <- c('vismodel', 'data.frame')
plot(colspace(eg1), col=rep(palette[1:2], each=50), 
     theta=155, phi=20, out.lwd=2, vert.cex=2, margin=c(0,0.5,0.5,0))
```

    ## Warning: Quantum catch are not relative, and have been transformed

``` r
text(x=grconvertX(-0.17,"npc"), y=grconvertY(1.15, "npc"), cex=1.5, "A") 

eg2 <- simulatedata[[1999]]
attr(eg2, 'conenumb') <- '4'
attr(eg2, 'relative') <- FALSE
attr(eg2, 'qcatch') <- 'Qi'
class(eg2) <- c('vismodel', 'data.frame')
plot(colspace(eg2), col=rep(palette[1:2], each=50), 
     theta=155, phi=20, out.lwd=2, vert.cex=2, margin=c(0,0.5,0.5,0))
```

    ## Warning: Quantum catch are not relative, and have been transformed

``` r
text(x=grconvertX(-0.17,"npc"), y=grconvertY(1.15, "npc"), cex=1.5, "B") 

dev.off()
```

    ## pdf 
    ##   2

``` r
######################
# RESULTS FROM SIM 1 #
######################
pdf(height=4*1.3, width=7*1.3, file='figures/samplesize_3.pdf')
par(mfrow=c(1,2), cex.lab=1.3, cex.axis=1.15, mar=c(5,4.5,4,1.5)+0.1)

plot(centdist~mahd, pch=21, 
     xlim=c(0.05, 10), ylim=c(0.01,10),
     col=NA, 
     bg=as.character(factor(adonisP, labels=palette[1:2])),
     log='xy', yaxt='n', xaxt='n',
     ylab='Mean distance (JND)', xlab='Effect size (Mahalanobis distance)')

axis(1, at=c(0.1, 1, 10), labels=c(0.1, 1, 10))
axis(1, at=c(seq(0.05,0.09, by=0.01), seq(0.2,0.9, by=0.1), seq(2,9, by=1)), tcl=par("tcl")*0.5, labels=FALSE)
axis(2, at=c(0.01,0.1, 1, 10), labels=c(0.01, 0.1, 1, 10))
axis(2, at=c(seq(0.02,0.09, by=0.01), seq(0.2,0.9, by=0.1), seq(2,9, by=1)), tcl=par("tcl")*0.5, labels=FALSE)

abline(h=1,lty=3, lwd=2)

#legend('topleft', 
legend(x=grconvertX(0,"npc"), y=grconvertY(0.93, "npc"),
       bty='n', pch=21, col=NA, pt.bg=palette[1:2], 
       legend=c('p > 0.05', 'p < 0.05'), cex=1.15)
text(x=grconvertX(0.05,"npc"), y=grconvertY(0.95, "npc"), cex=1.5, "A") 

plot(overlap~mahd, pch=21, 
     xlim=c(0.05, 10), ylim=c(0,80),
     bg=as.character(factor(adonisP, labels=palette[1:2])), 
     col=NA,
     log='x', yaxt='n', xaxt='n',
     ylab='Colour volume overlap (%)', xlab='Effect size (Mahalanobis distance)')

axis(1, at=c(0.1, 1, 10), labels=c(0.1, 1, 10))
axis(1, at=c(seq(0.06,0.09, by=0.01), seq(0.2,0.9, by=0.1), seq(2,9, by=1)), tcl=par("tcl")*0.5, labels=FALSE)
axis(2, at=c(0, 20, 40, 60, 80))

text(x=grconvertX(0.05,"npc"), y=grconvertY(0.95, "npc"), cex=1.5, "B") 

#par(fig=c(0.78,1, 0,1), new = TRUE)
#plot(as.numeric(tapply(adonisP, cut(overlap, seq(0,100, by=10)), mean)), 
#     I(seq(0,100, by=10) + 5)[-11], 
#     pch=21, col=brewer.pal(5,'PRGn')[5], bg=brewer.pal(5,'PRGn')[5], 
#     type='b', xlim=c(-0.1,1.1), ylim=c(0,80), cex=1.5, lwd=1.5, 
#     xlab="", xaxt='n', ylab="", yaxt='n')
#axis(1, at=c(0,1))
#axis(1, at=c(0.2,0.4,0.6,0.8), tcl=par("tcl")*0.5, labels=FALSE)

dev.off()
```

    ## pdf 
    ##   2
