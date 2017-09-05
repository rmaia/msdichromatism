Worked examples and different approaches
================

-   [For a simple example looking at a single patch, [click here](lizardexample-single.md)](#for-a-simple-example-looking-at-a-single-patch-click-here)
-   [The Data](#the-data)
-   [Approach 1: distance-based PERMANOVA](#approach-1-distance-based-permanova)
    -   [Testing for separation among groups](#testing-for-separation-among-groups)
    -   [Testing for above-threshold mean differences between groups](#testing-for-above-threshold-mean-differences-between-groups)
-   [Approach 2: Cartesian MANOVA](#approach-2-cartesian-manova)
    -   [Testing for separation among groups](#testing-for-separation-among-groups-1)
    -   [Testing for above-threshold mean differences between groups](#testing-for-above-threshold-mean-differences-between-groups-1)
-   [Approach 3: Bayesian multi-response model using MCMCglmm](#approach-3-bayesian-multi-response-model-using-mcmcglmm)
    -   [Testing for separation among groups](#testing-for-separation-among-groups-2)
    -   [Testing for above-threshold mean differences between groups](#testing-for-above-threshold-mean-differences-between-groups-2)
-   [Comparing Bayesian vs. Bootstrap mean difference estimates](#comparing-bayesian-vs.-bootstrap-mean-difference-estimates)

#### For a simple example looking at a single patch, [click here](lizardexample-single.md)

First, we need to install the bleeding edge version of pavo:

``` r
devtools::install_github('rmaia/pavo@jnd2xyz')
```

load the necessary packages and functions:

``` r
require(pavo)
require(MCMCglmm)
require(scatterplot3d)
require(gridExtra)
require(vegan)
require(MASS)
require(RColorBrewer)

# load aesthetic functions (plot, make colors transparent)
source('R/aesthetic.R')

# load function to convert JND to cartesian coordinates
source('R/jnd2xyz.R')

# load function for bootstrap
source('R/bootstrapcentroiddS.R')

# Distance matrix generator
distmat <- function(x){
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
  M
  }

#color palette
palette <- rcbalpha(0.8, 4, 'Set1')
```

The Data
========

Reflectance data from four body regions of male and female *Ctenophorus ornatus* (Whiting et al. 2015, Biol J Linn Soc). Labium, throat, tongue, and mouth-roof.

**Question:** Which body regions are sexually dichromatic?

Calculate deltaS according to conspecific (tetrachromatic) visual system

``` r
specs <- list(
              lab = as.rspec(read.csv('data/dichromatism/lab.csv'), interp = FALSE),
              throat = as.rspec(read.csv('data/dichromatism/throat.csv'), interp = FALSE),
              roof = as.rspec(read.csv('data/dichromatism/roof.csv'), interp = FALSE),
              tongue = as.rspec(read.csv('data/dichromatism/tongue.csv'), interp = FALSE)
              )

# Ctenophorus ornatus
liz_vis <- sensmodel(c(360, 440, 493, 571)) 
names(liz_vis) <- c('wl', 'u', 's', 'm', 'l')

models <- lapply(specs, vismodel, visual = liz_vis, relative = FALSE, qcatch='Qi')

spaces <- lapply(models, colspace)

deltaS <- lapply(models, coldist, achro = FALSE, n = c(1,1,3.5,6), 
                                  weber = 0.1, noise = "neural")
```

Visualise

``` r
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 4, 2, byrow = TRUE))

aggplot(specs[['lab']], by=gsub("[0-9].*","",names(specs[['lab']])), lwd=3, ylim=c(0,50))
text(x=grconvertX(0.15,"npc"), y=grconvertY(0.8, "npc"), cex=1.5, pos=3, "Labials") 

scatterplot3d(spaces[['lab']][,c('x','y','z')],  
      bg=as.character(factor(gsub("[0-9].*","",names(specs[['lab']]))[-1], 
                             labels=palette[1:2])),
      box=FALSE, pch=21, cex.symbols=2, color=NA,
      x.ticklabs='', y.ticklabs='', z.ticklabs='', xlab='', ylab='', zlab='')


aggplot(specs[['throat']], by=gsub("[0-9].*","",names(specs[['throat']])), lwd=3, ylim=c(0,50))
text(x=grconvertX(0.15,"npc"), y=grconvertY(0.8, "npc"), cex=1.5, pos=3, "Throat")

scatterplot3d(spaces[['throat']][,c('x','y','z')],  
      bg=as.character(factor(gsub("[0-9].*","",names(specs[['throat']]))[-1],
                             labels=palette[1:2])),
      box=FALSE, pch=21, cex.symbols=2, color=NA,
      x.ticklabs='', y.ticklabs='', z.ticklabs='', xlab='', ylab='', zlab='')


aggplot(specs[['roof']], by=gsub("[0-9].*","",names(specs[['roof']])), lwd=3, ylim=c(0,50))
text(x=grconvertX(0.15,"npc"), y=grconvertY(0.8, "npc"), cex=1.5, pos=3, "Roof")

scatterplot3d(spaces[['roof']][,c('x','y','z')],  
      bg=as.character(factor(gsub("[0-9].*","",names(specs[['roof']]))[-1],
                             labels=palette[1:2])),
      box=FALSE, pch=21, cex.symbols=2, color=NA,
      x.ticklabs='', y.ticklabs='', z.ticklabs='', xlab='', ylab='', zlab='')


aggplot(specs[['tongue']], by=gsub("[0-9].*","",names(specs[['tongue']])), lwd=3, ylim=c(0,50))
text(x=grconvertX(0.15,"npc"), y=grconvertY(0.8, "npc"), cex=1.5, pos=3, "Tongue")

scatterplot3d(spaces[['tongue']][,c('x','y','z')],  
      bg=as.character(factor(gsub("[0-9].*","",names(specs[['tongue']]))[-1],
                             labels=palette[1:2])),
      box=FALSE, pch=21, cex.symbols=2, color=NA,
      x.ticklabs='', y.ticklabs='', z.ticklabs='', xlab='', ylab='', zlab='')
```

![](../output/figures/examples/lizard_fig_liz_spectcs-1.jpeg)

Approach 1: distance-based PERMANOVA
====================================

``` r
# Setup distance matrices & groupings for each body part
mat <- list(
            lab = distmat(deltaS$lab),
            throat = distmat(deltaS$throat),
            roof = distmat(deltaS$roof),
            tongue = distmat(deltaS$tongue)
            )

group <- lapply(mat, function(x) substring(rownames(as.matrix(x)), 1, 1))
```

Testing for separation among groups
-----------------------------------

Fist, let's test the assumption of homogeneity of variances

``` r
bdisp <- lapply(names(mat), function(x) betadisper(mat[[x]], group[[x]]))
names(bdisp) <- names(mat)
lapply(bdisp, anova)
```

    ## $lab
    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value   Pr(>F)   
    ## Groups     1  2.7904 2.79044  8.1507 0.005993 **
    ## Residuals 57 19.5143 0.34236                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $throat
    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Groups     1  0.031 0.03125   0.037 0.8482
    ## Residuals 58 49.040 0.84552               
    ## 
    ## $roof
    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Groups     1  1.363 1.36296  2.1709 0.1463
    ## Residuals 55 34.531 0.62784               
    ## 
    ## $tongue
    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## Groups     1    3.1 3.10022  5.0227 0.02886 *
    ## Residuals 58   35.8 0.61724                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Labials and tongue have unequal variances, which can influence the results of the PERMANOVA. When this assumption is violated but the largest group has the largest variance, it usually isn't a big problem. That is the case for tongue...

``` r
TukeyHSD(bdisp$tongue)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##          diff        lwr       upr     p adj
    ## M-F 0.4569122 0.04881091 0.8650135 0.0288602

``` r
# Sample sizes
table(group$tongue)
```

    ## 
    ##  F  M 
    ## 27 33

...but not for labials:

``` r
TukeyHSD(bdisp$lab)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##           diff        lwr        upr     p adj
    ## M-F -0.4365211 -0.7426987 -0.1303434 0.0059932

``` r
# Sample sizes
table(group$lab)
```

    ## 
    ##  F  M 
    ## 27 32

So the latter's PERMANOVA should be treated with caution.

Permutational MANOVA:

``` r
pmanova <- lapply(names(mat), function(x) adonis(mat[[x]] ~ group[[x]]))
names(pmanova) <- names(mat)
pmanova
```

    ## $lab
    ## 
    ## Call:
    ## adonis(formula = mat[[x]] ~ group[[x]]) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## group[[x]]  1    24.509 24.5094  13.964 0.19678  0.001 ***
    ## Residuals  57   100.042  1.7551         0.80322           
    ## Total      58   124.552                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $throat
    ## 
    ## Call:
    ## adonis(formula = mat[[x]] ~ group[[x]]) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## group[[x]]  1    42.565  42.565  14.842 0.20376  0.001 ***
    ## Residuals  58   166.335   2.868         0.79624           
    ## Total      59   208.901                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $roof
    ## 
    ## Call:
    ## adonis(formula = mat[[x]] ~ group[[x]]) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## group[[x]]  1     0.563 0.56278 0.52282 0.00942  0.495
    ## Residuals  55    59.203 1.07642         0.99058       
    ## Total      56    59.766                 1.00000       
    ## 
    ## $tongue
    ## 
    ## Call:
    ## adonis(formula = mat[[x]] ~ group[[x]]) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## group[[x]]  1     2.162  2.1617  1.6293 0.02732  0.223
    ## Residuals  58    76.953  1.3268         0.97268       
    ## Total      59    79.115                 1.00000

Testing for above-threshold mean differences between groups
-----------------------------------------------------------

We will use a bootstrap approach to consider the uncertainty on the mean difference between groups when comparing it to a threshold value of 1 JND. Normally it would be unnecessary to do this for Mouth-Roof and Tongue, since the test of separation indicates these are not discriminable. Therefore, even if their mean differences are above 1 JND, males and females would not be considered different, and it would just indicate that within-group variation is also high. But for illustrative purposes, we will analyze all body parts.

``` r
# Groups
models$lab$group <- substring(rownames(models$lab), 1, 1)
models$throat$group <- substring(rownames(models$throat), 1, 1)
models$roof$group <- substring(rownames(models$roof), 1, 1)
models$tongue$group <- substring(rownames(models$tongue), 1, 1)

bootds <- lapply(models, function(x) bootcentroidDS(x[,1:4], x[,5], n=c(1,1,3.5,6), weber=0.1, qcatch='Qi', achro=FALSE))
```

    ## Warning: number of cones not specified; assumed to be 4

    ## Warning: number of cones not specified; assumed to be 4

    ## Warning: number of cones not specified; assumed to be 4

    ## Warning: number of cones not specified; assumed to be 4

``` r
bootds
```

    ## $lab
    ##     measured.dS    CI.lwr   CI.upr
    ## F-M    1.293707 0.8727416 1.743769
    ## 
    ## $throat
    ##     measured.dS   CI.lwr   CI.upr
    ## F-M    1.693031 1.119148 2.327343
    ## 
    ## $roof
    ##     measured.dS     CI.lwr    CI.upr
    ## F-M   0.1990052 0.04603181 0.7221358
    ## 
    ## $tongue
    ##     measured.dS    CI.lwr    CI.upr
    ## F-M   0.3815314 0.1280804 0.8118222

We can see that, though labium is statistically significant, the distance between groups cannot be considered to be above threshold:

``` r
bootres <- do.call(rbind, bootds) 
rownames(bootres) <- c('Labium', 'Throat', 'Roof', 'Tongue')

plot(bootres[,1], xlim=c(0.5, 4.5), ylim=c(0, 2.5), pch=21, bg=1, cex=2, xaxt='n', xlab='Centroid comparison', ylab='Chromatic contrast (JND)')

abline(h=1, lty=3, lwd=2)
segments(1:4, bootres[,2], 1:4, bootres[,3], lwd=2)

axis(1, at=1:4, labels=rownames(bootres))
```

![](../output/figures/examples/lizard_fig_unnamed-chunk-7-1.jpeg)

Approach 2: Cartesian MANOVA
============================

First, we need to convert our distances into perceptually-proportional cartesian coordinates:

``` r
# get perceptual xyz and group vector
pxyz <- lapply(names(deltaS), function(x) data.frame(
  jnd2xyz(deltaS[[x]]),                           # perceptual xyz
  group = substring(names(specs[[x]]), 1, 1)[-1]) # group vector
  )

names(pxyz) <- names(deltaS)
```

Testing for separation among groups
-----------------------------------

... add test of variances ...

run MANOVA on cartesian coordinates

``` r
lapply(pxyz, function(x) summary(manova(lm(cbind(x,y,z)~group, data = x))))
```

    ## $lab
    ##           Df  Pillai approx F num Df den Df    Pr(>F)    
    ## group      1 0.43847   14.316      3     55 5.201e-07 ***
    ## Residuals 57                                             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $throat
    ##           Df Pillai approx F num Df den Df    Pr(>F)    
    ## group      1 0.4017   12.533      3     56 2.229e-06 ***
    ## Residuals 58                                            
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $roof
    ##           Df   Pillai approx F num Df den Df Pr(>F)
    ## group      1 0.049607  0.92213      3     53 0.4365
    ## Residuals 55                                       
    ## 
    ## $tongue
    ##           Df   Pillai approx F num Df den Df Pr(>F)
    ## group      1 0.077482   1.5678      3     56 0.2073
    ## Residuals 58

Testing for above-threshold mean differences between groups
-----------------------------------------------------------

We will use a bootstrap approach to consider the uncertainty on the mean difference between groups when comparing it to a threshold value of 1 JND. Normally it would be unnecessary to do this for Mouth-Roof and Tongue, since the test of separation indicates these are not discriminable. Therefore, even if their mean differences are above 1 JND, males and females would not be considered different, and it would just indicate that within-group variation is also high. But for illustrative purposes, we will analyze all body parts.

``` r
# Groups
models$lab$group <- substring(rownames(models$lab), 1, 1)
models$throat$group <- substring(rownames(models$throat), 1, 1)
models$roof$group <- substring(rownames(models$roof), 1, 1)
models$tongue$group <- substring(rownames(models$tongue), 1, 1)

bootds <- lapply(models, function(x) bootcentroidDS(x[,1:4], x[,5], n=c(1,1,3.5,6), weber=0.1, qcatch='Qi', achro=FALSE))
```

    ## Warning: number of cones not specified; assumed to be 4

    ## Warning: number of cones not specified; assumed to be 4

    ## Warning: number of cones not specified; assumed to be 4

    ## Warning: number of cones not specified; assumed to be 4

``` r
bootds
```

    ## $lab
    ##     measured.dS    CI.lwr   CI.upr
    ## F-M    1.293707 0.8805381 1.751385
    ## 
    ## $throat
    ##     measured.dS   CI.lwr   CI.upr
    ## F-M    1.693031 1.142767 2.363832
    ## 
    ## $roof
    ##     measured.dS    CI.lwr    CI.upr
    ## F-M   0.1990052 0.0558494 0.7005666
    ## 
    ## $tongue
    ##     measured.dS    CI.lwr    CI.upr
    ## F-M   0.3815314 0.1174271 0.8522395

We can see that, though labium is statistically significant, the distance between groups cannot be considered to be above threshold:

``` r
bootres <- do.call(rbind, bootds) 
rownames(bootres) <- c('Labium', 'Throat', 'Roof', 'Tongue')

plot(bootres[,1], xlim=c(0.5, 4.5), ylim=c(0, 2.5), pch=21, bg=1, cex=2, xaxt='n', xlab='Centroid comparison', ylab='Chromatic contrast (JND)')

abline(h=1, lty=3, lwd=2)
segments(1:4, bootres[,2], 1:4, bootres[,3], lwd=2)

axis(1, at=1:4, labels=rownames(bootres))
```

![](../output/figures/examples/lizard_fig_unnamed-chunk-10-1.jpeg)

Approach 3: Bayesian multi-response model using MCMCglmm
========================================================

we can extract both statistical and perceptual information from the posterior distribution of a bayesian analysis using a multi-response model and the package MCMCglmm. In this analysis, we will estimate the variance-covariance structure, so the assumption of homogeneity of variances is relaxed.

Before we do that, we will center all cartesian variables on the female means --- this way, since the model does not have an intercept, male estimates and pMCMC can be interpreted as the effect of sex on each response (X, Y and Z). **This step is essential if we want to interpret model estimates as difference estimates**:

``` r
# get perceptual xyz and group vector
pxyz <- lapply(names(deltaS), function(x) data.frame(
  jnd2xyz(deltaS[[x]]),                           # perceptual xyz
  group = substring(names(specs[[x]]), 1, 1)[-1]) # group vector
  )

names(pxyz) <- names(deltaS)

# Centering variables to the mean of the first group (in this case, F)
pxyz$lab[, -4] <- sweep(as.matrix(pxyz$lab[,-4]), 2, 
                        as.matrix(aggregate(pxyz$lab[,-4], list(pxyz$lab[,4]), mean)[1,-1]), '-') 

pxyz$throat[, -4] <- sweep(as.matrix(pxyz$throat[,-4]), 2, 
                        as.matrix(aggregate(pxyz$throat[,-4], list(pxyz$throat[,4]), mean)[1,-1]), '-') 

pxyz$roof[, -4] <- sweep(as.matrix(pxyz$roof[,-4]), 2, 
                        as.matrix(aggregate(pxyz$roof[,-4], list(pxyz$roof[,4]), mean)[1,-1]), '-') 

pxyz$tongue[, -4] <- sweep(as.matrix(pxyz$tongue[,-4]), 2, 
                        as.matrix(aggregate(pxyz$tongue[,-4], list(pxyz$tongue[,4]), mean)[1,-1]), '-') 
```

Testing for separation among groups
-----------------------------------

Below we will run the models, but for simplicity we will only show the verbose output for the tests on the labials.

``` r
# Running MCMCglmm models with default priors
mcmcres <- lapply(pxyz, function(x) MCMCglmm(
  fixed = cbind(x,y,z) ~ trait:group - 1,
  rcov = ~us(trait):units,
  family = rep('gaussian', 3),
  data = x,
  nit=11000, burnin=1000, thin=10,
  verbose = FALSE
)
)
```

Chains seem to be mixing ok-ish:

``` r
plot(mcmcres[['lab']]$VCV, density = FALSE)
```

![](../output/figures/examples/lizard_fig_unnamed-chunk-12-1.jpeg)

``` r
plot(mcmcres[['lab']]$Sol, density = FALSE)
```

![](../output/figures/examples/lizard_fig_unnamed-chunk-12-2.jpeg)

We can see that the posterior estimates approximate the data fairly well (though the chain could probably have been run for a little longer):

``` r
summary(mcmcres[['lab']])
```

    ## 
    ##  Iterations = 1001:10991
    ##  Thinning interval  = 10
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 284.307 
    ## 
    ##  R-structure:  ~us(trait):units
    ## 
    ##                     post.mean  l-95% CI u-95% CI eff.samp
    ## traitx:traitx.units   0.70071  0.452601   0.9624     1000
    ## traity:traitx.units  -0.03001 -0.258707   0.2227     1122
    ## traitz:traitx.units   0.10702  0.006799   0.2132     1000
    ## traitx:traity.units  -0.03001 -0.258707   0.2227     1122
    ## traity:traity.units   1.01992  0.680660   1.4255     1000
    ## traitz:traity.units  -0.36922 -0.522232  -0.2200     1000
    ## traitx:traitz.units   0.10702  0.006799   0.2132     1000
    ## traity:traitz.units  -0.36922 -0.522232  -0.2200     1000
    ## traitz:traitz.units   0.17907  0.114523   0.2494     1000
    ## 
    ##  Location effects: cbind(x, y, z) ~ trait:group - 1 
    ## 
    ##                post.mean   l-95% CI   u-95% CI eff.samp  pMCMC    
    ## traitx:groupF  0.0004106 -0.3043810  0.3103640     1000  0.994    
    ## traity:groupF -0.0012237 -0.3804630  0.3744444     1000  1.000    
    ## traitz:groupF  0.0001009 -0.1518630  0.1807687     1000  0.984    
    ## traitx:groupM -1.2473926 -1.5325616 -0.9307164     1057 <0.001 ***
    ## traity:groupM -0.3351477 -0.7010905  0.0304231     1000  0.072 .  
    ## traitz:groupM  0.0687383 -0.0882124  0.2150931     1000  0.392    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# checking if model estimates approximate empirical means
# First column = Females, second column = Males
# (remember: we centered on Female means, so all female estimates should be zero)
round(as.matrix(aggregate(pxyz$lab[,c('x','y','z')], by= list(pxyz$lab[,'group']), mean)[,-1]), 4)
```

    ##           x       y      z
    ## [1,]  0.000  0.0000 0.0000
    ## [2,] -1.249 -0.3306 0.0662

``` r
round(matrix(summary(mcmcres[['lab']]$Sol)$statistics[,'Mean'], nrow=2, byrow=TRUE), 4)
```

    ##         [,1]    [,2]   [,3]
    ## [1,]  0.0004 -0.0012 0.0001
    ## [2,] -1.2474 -0.3351 0.0687

The results above already indicate that the groups are statistically different along one of the axes (x). Some people might consider this to be subject to issues of multiple comparisons, because we are considering groups to be different if they are different in **any** (i.e. at least 1) dimension, not **all** dimensions together.

We don't agree, because the distributions of each response variables are not being calculated independently. Their multivariate nature of the responses is accounted for and the model is estimating their covariances --- so it's **not** the same as, say, testing each of them individually.

Nonetheless, we can obtain an empirical P-value from the MCMC chain that is analogous to a Wald test, that is, testing if **all** responses are **jointly** different than zero, by estimating the ellipsoid representing its joint credible interval and seeing how many iterations of the MCMC chain include zero:

``` r
# function modified from: https://stat.ethz.ch/pipermail/r-help/2006-September/113184.html

multimcmcpval <- function(samp)
{
   ## elementary version that creates an empirical p-value for the
   ## hypothesis that the columns of samp have mean zero versus a
   ## general multivariate distribution with elliptical contours.

   ## differences from the mean standardized by the observed
   ## variance-covariance factor
   std <- backsolve(chol(var(samp)),
                    cbind(0, t(samp)) - colMeans(samp),
                    transpose = TRUE)
   sqdist <- colSums(std * std)
   #sum(sqdist[-1] > sqdist[1])/nrow(samp)
   mean(sqdist[-1] > sqdist[1])
}

lapply(mcmcres, function(x) multimcmcpval(x$Sol[,4:6]))
```

    ## $lab
    ## [1] 0
    ## 
    ## $throat
    ## [1] 0
    ## 
    ## $roof
    ## [1] 0.127
    ## 
    ## $tongue
    ## [1] 0.028

Interestingly, the tongue is coming out as significantly different, likely a result of the differences in distribution between males and females (see first figure).

Testing for above-threshold mean differences between groups
-----------------------------------------------------------

We will use the posterior distributions of estimates from male and female xyz values to estimate their mean differences along with the credible intervals along these estimates, comparing them to a threshold value of 1 JND.

``` r
dmcmc <- lapply(mcmcres, function(y) apply(
  y$Sol, 1, function(x) as.numeric(dist(rbind(x[1:3], x[4:6])))
  )
  )


credibleints <- cbind(
do.call(rbind, lapply(dmcmc, function(x) posterior.mode(mcmc(x)))),
do.call(rbind, lapply(dmcmc, function(x) HPDinterval(mcmc(x))))
)

colnames(credibleints)[1] <- 'posteriormode'

credibleints
```

    ##        posteriormode      lower     upper
    ## lab        1.3399752 0.94894244 1.7564114
    ## throat     1.7572383 1.11434816 2.4488245
    ## roof       0.1635817 0.02179937 0.6988794
    ## tongue     0.2985872 0.09552225 0.8230676

``` r
par(mfrow=c(2,2))
dplots <- lapply(names(dmcmc), function(x) {plot(density(dmcmc[[x]]), main=x); abline(v=1, lty=3)})
```

![](../output/figures/examples/lizard_fig_unnamed-chunk-15-1.jpeg)

Comparing Bayesian vs. Bootstrap mean difference estimates
==========================================================

We can compare what we obtain from the Bayesian approach with the results from the Bootstrap:

``` r
palette <- rcbalpha(1, 4, 'Set1')

plot(y=bootres[,'measured.dS'], x=c(0.9, 1.9, 2.9, 3.9), xlim=c(0.5,4.5), ylim=c(0, 2.5), pch=21, col=NULL, bg=palette[3], cex=2, xaxt='n', xlab='Centroid comparison', ylab='Chromatic contrast (JND)')

segments(c(0.9, 1.9, 2.9, 3.9), bootres[,'CI.lwr'], c(0.9, 1.9, 2.9, 3.9), bootres[,'CI.upr'], lwd=2, col=palette[3])

points(y=credibleints[,'posteriormode'], x=c(1.1, 2.1, 3.1, 4.1), pch=21, bg=palette[4], cex=2, col=NULL)

segments(c(1.1, 2.1, 3.1, 4.1), credibleints[,'lower'], c(1.1, 2.1, 3.1, 4.1), credibleints[,'upper'], lwd=2, col=palette[4])

abline(h=1, lty=3, lwd=2)


axis(1, at=1:4, labels=rownames(bootres))

legend('topright', pch=21, pt.cex=1.6, col=palette[3:4], lwd=2, pt.bg=palette[3:4], legend=c('Bootsrap estimate', 'Bayesian posterior'))
```

![](../output/figures/examples/lizard_fig_unnamed-chunk-16-1.jpeg)
