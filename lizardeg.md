Example w/ real data
====================

### Reflectance data from several body regions of male and female *Ctenophorus ornatus* (Whiting et al. 2015, Biol J Linn Soc)

Calculate deltaS

    ## wavelengths found in column 1 
    ## wavelengths found in column 1 
    ## wavelengths found in column 1 
    ## wavelengths found in column 1

Plot 'em

``` r
par(pty="s", mfrow = c(2, 2))

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$lab[grepl("M", rownames(models_rel$lab)), ])[, c('x','y','z')]), pch=19,
                      xlim=c(-0.05, 0.05), ylim=c(-0.05,0.02), zlim=c(-0.1, 0.4), box=F, main = 'labium')
sp3d$points3d(suppressWarnings(tcs(models_rel$lab[grepl("F", rownames(models_rel$lab)), ])[, c('x','y','z')]), col='red',pch=19)

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$throat[grepl("M", rownames(models_rel$throat)), ])[, c('x','y','z')]), pch=19,
                      xlim=c(-0.05, 0.05), ylim=c(-0.05,0.02), zlim=c(-0.1, 0.4), box=F, main = 'throat')
sp3d$points3d(suppressWarnings(tcs(models_rel$throat[grepl("F", rownames(models_rel$throat)), ])[, c('x','y','z')]), col='red',pch=19)

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$roof[grepl("M", rownames(models_rel$roof)), ])[, c('x','y','z')]), pch=19,
                      xlim=c(-0.05, 0.05), ylim=c(-0.05,0.02), zlim=c(-0.1, 0.4), box=F, main = 'roof')
sp3d$points3d(suppressWarnings(tcs(models_rel$roof[grepl("F", rownames(models_rel$roof)), ])[, c('x','y','z')]), col='red',pch=19)

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$tongue[grepl("M", rownames(models_rel$tongue)), ])[, c('x','y','z')]), pch=19,
                      xlim=c(-0.05, 0.05), ylim=c(-0.05,0.02), zlim=c(-0.1, 0.4), box=F, main = 'tongue')
sp3d$points3d(suppressWarnings(tcs(models_rel$tongue[grepl("F", rownames(models_rel$tongue)), ])[, c('x','y','z')]), col='red',pch=19)
```

![](output/figures/lizardeg/lizardeg_figtcs-1.png)<!-- -->

``` r
# Check 'em out
p1 <- ggplot(deltaS$lab, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
        facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
        ggtitle('labial') + theme(legend.position="none")

p2 <- ggplot(deltaS$throat, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
        facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
        ggtitle('throat') + theme(legend.position="none")

p3 <- ggplot(deltaS$roof, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
        facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
        ggtitle('roof') + theme(legend.position="none")

p4 <- ggplot(deltaS$tongue, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
        facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
        ggtitle('tongue') + theme(legend.position="none")

grid.arrange(p1, p2, p3, p4, ncol=2)
```

![](output/figures/lizardeg/lizardeg_figdeltaplot-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 3.2.3 (2015-12-10)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: OS X 10.11.5 (El Capitan)
    ## 
    ## locale:
    ## [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] gridExtra_2.0.0      ggplot2_2.0.0        scatterplot3d_0.3-37
    ## [4] pavo_0.5-5           rgl_0.95.1441       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.3      knitr_1.12.3     magrittr_1.5     maps_3.1.0      
    ##  [5] magic_1.5-6      munsell_0.4.3    colorspace_1.2-6 geometry_0.3-6  
    ##  [9] plyr_1.8.3       stringr_1.0.0    tools_3.2.3      grid_3.2.3      
    ## [13] gtable_0.1.2     htmltools_0.3    yaml_2.1.13      digest_0.6.9    
    ## [17] reshape2_1.4.1   mapproj_1.2-4    formatR_1.2.1    rcdd_1.1-10     
    ## [21] evaluate_0.8     rmarkdown_0.9.5  labeling_0.3     stringi_1.0-1   
    ## [25] scales_0.3.0
