(can someone tell me how the hell to choose the folder to put the files resulting from building this markdown? UGH)

Simulation
==========

Here we'll attempt to show, using simulations, that using mean deltaS to estimate the difference between two groups might not be appropriate.

``` r
require(pavo)
```

    ## Loading required package: pavo

    ## Loading required package: rgl

``` r
require(scatterplot3d)
```

    ## Loading required package: scatterplot3d

``` r
require(ggplot2)
```

    ## Loading required package: ggplot2

``` r
sessionInfo()
```

    ## R version 3.3.0 (2016-05-03)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: OS X 10.11.4 (El Capitan)
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_2.1.0        scatterplot3d_0.3-37 pavo_0.5-5          
    ## [4] rgl_0.95.1441       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.5        knitr_1.13         magrittr_1.5      
    ##  [4] maps_3.1.0         magic_1.5-6        munsell_0.4.3     
    ##  [7] colorspace_1.2-6   geometry_0.3-6     plyr_1.8.3        
    ## [10] stringr_1.0.0      tools_3.3.0        grid_3.3.0        
    ## [13] gtable_0.2.0       htmltools_0.3.5    yaml_2.1.13       
    ## [16] digest_0.6.9       mapproj_1.2-4      rcdd_1.1-10       
    ## [19] evaluate_0.9       rmarkdown_0.9.6.10 stringi_1.0-1     
    ## [22] scales_0.4.0

### Example 1: low intra-group variability, low inter-group distance

``` r
# step 1: generate data

set.seed(1982)

# we'll just consider usml are uncorrelated for this example

# we'll generate data from a lognormal distribution to avoid negative values
# and a variance (on the log scale) of 0.002

groupA <- data.frame(
  u = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.002)),
  s = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.002)),
  m = rlnorm(30, meanlog=log(0.3), sdlog=sqrt(0.002)),
  l = rlnorm(30, meanlog=log(0.7), sdlog=sqrt(0.002))
)

groupB <- data.frame(
  u = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.002)),
  s = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.002)),
  m = rlnorm(30, meanlog=log(0.3), sdlog=sqrt(0.002)),
  l = rlnorm(30, meanlog=log(0.7), sdlog=sqrt(0.002))
)


colnames(groupA) <- colnames(groupB) <- c('u','s','m', 'l')
attr(groupA, 'relative') <- attr(groupB, 'relative') <- FALSE

par(pty="s")
sp3d <- scatterplot3d(suppressWarnings(tcs(groupA)[, c('x','y','z')]), pch=19,
                      xlim=c(0.26,0.38), ylim=c(-0.17,0.024), zlim=c(-0.2,-0.14), box=F)
sp3d$points3d(suppressWarnings(tcs(groupB)[, c('x','y','z')]), col='red',pch=19)
```

![](Simulations_files/figure-markdown_github/unnamed-chunk-2-1.png)

Note that although USML were simulated uncorrelated, XYZ are correlated.

Calculate deltaS

``` r
rownames(groupA) <- paste('gA',1:30,sep='')
rownames(groupB) <- paste('gB',1:30,sep='')

alldat <- rbind(groupA, groupB)

deltaS <- coldist(alldat, achro=FALSE)

deltaS$comparison <- NA

deltaS$comparison[grepl('A', deltaS$patch1) & grepl('A', deltaS$patch2)] <- 'intra.A'
deltaS$comparison[grepl('B', deltaS$patch1) & grepl('B', deltaS$patch2)] <- 'intra.B'
deltaS$comparison[grepl('A', deltaS$patch1) & grepl('B', deltaS$patch2)] <- 'inter'

ggplot(deltaS, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
  facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
  theme(legend.position="none")
```

![](Simulations_files/figure-markdown_github/unnamed-chunk-3-1.png)

### Example 2: High intra-group variability, low inter-group distance

``` r
# step 1: generate data

set.seed(1982)

# we'll just consider usml are uncorrelated for this example

# we'll generate data from a lognormal distribution to avoid negative values
# and a variance (on the log scale) of 0.002

groupA <- data.frame(
  u = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.02)),
  s = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.02)),
  m = rlnorm(30, meanlog=log(0.3), sdlog=sqrt(0.02)),
  l = rlnorm(30, meanlog=log(0.7), sdlog=sqrt(0.02))
)

groupB <- data.frame(
  u = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.02)),
  s = rlnorm(30, meanlog=log(0.1), sdlog=sqrt(0.02)),
  m = rlnorm(30, meanlog=log(0.3), sdlog=sqrt(0.02)),
  l = rlnorm(30, meanlog=log(0.7), sdlog=sqrt(0.02))
)


colnames(groupA) <- colnames(groupB) <- c('u','s','m', 'l')
attr(groupA, 'relative') <- attr(groupB, 'relative') <- FALSE

par(pty="s")
sp3d <- scatterplot3d(suppressWarnings(tcs(groupA)[, c('x','y','z')]), pch=19,
                      xlim=c(0.26,0.38), ylim=c(-0.17,0.024), zlim=c(-0.2,-0.14), box=F)
sp3d$points3d(suppressWarnings(tcs(groupB)[, c('x','y','z')]), col='red',pch=19)
```

![](Simulations_files/figure-markdown_github/unnamed-chunk-4-1.png)

Note that although USML were simulated uncorrelated, XYZ are correlated.

Calculate deltaS

``` r
rownames(groupA) <- paste('gA',1:30,sep='')
rownames(groupB) <- paste('gB',1:30,sep='')

alldat <- rbind(groupA, groupB)

deltaS <- coldist(alldat, achro=FALSE)

deltaS$comparison <- NA

deltaS$comparison[grepl('A', deltaS$patch1) & grepl('A', deltaS$patch2)] <- 'intra.A'
deltaS$comparison[grepl('B', deltaS$patch1) & grepl('B', deltaS$patch2)] <- 'intra.B'
deltaS$comparison[grepl('A', deltaS$patch1) & grepl('B', deltaS$patch2)] <- 'inter'

ggplot(deltaS, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
  facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
  theme(legend.position="none")
```

![](Simulations_files/figure-markdown_github/unnamed-chunk-5-1.png)

#### Note that in both cases there is no dichromatism, but regular analyses would consider the second case dichromatic.

How to test this?
-----------------

I think something from the likes of the adonis function in vegan, but I don't know how to use a custom distance matrix with it yet.
