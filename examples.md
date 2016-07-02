Possible examples w/ real data
==============================

``` r
# Stolen from simspt2
adoniscoldist <- function(x){
  dmat <- matrix(0, nrow=length(unique(x$patch1)), ncol=length(unique(x$patch1)))
  rownames(dmat) <- colnames(dmat) <- as.character(unique(x$patch1))
  
  for(i in rownames(dmat))
    for(j in colnames(dmat))
      if(length(x$dS[x$patch1 == i & x$patch2 == j]) != 0)
      dmat[i,j] <- dmat[j,i] <- x$dS[x$patch1 == i & x$patch2 == j]
  
  grouping <- gsub('[0-9]','', rownames(dmat))
  
  adonis(dmat~grouping)
}

# Trichromat space & plot (from pavo)
trispace <- function(vismodeldata){
  
  dat <- vismodeldata
    
  s <- dat[, 1]
  m <- dat[, 2]
  l <- dat[, 3]
  
# cartesian coordinates
    
  x <- (1/sqrt(2)) * (l - m) 
  y <- (sqrt(2)/sqrt(3)) * (s - ((l + m)/2))   
    
# colorimetrics
  r.vec <- sqrt((abs(x))^2 + (abs(y)^2))
  h.theta <- atan2(y, x)
  
  res.p <- data.frame(s, m, l, x, y, h.theta, r.vec, row.names = rownames(dat))
  
  res <- res.p
  
  res
}

triplot <- function(tridata, labels = TRUE, achro = TRUE, achrocol = 'grey', achrosize = 0.8, 
                     cex.labels = 1, out.lwd = 1, out.lcol = 'black', out.lty = 1, ...){ 
  
  arg <- list(...)
  
# Set defaults
  if(is.null(arg$col))
    arg$col <- 'black'
  if(is.null(arg$pch))
    arg$pch <- 19
  if(is.null(arg$type))
    arg$type = 'p'
  if(is.null(arg$xlim))
    arg$xlim <- c(-1/sqrt(2), 1/sqrt(2))
  if(is.null(arg$ylim))
    arg$ylim <- c(-sqrt(2)/(2*(sqrt(3))), sqrt(2)/sqrt(3))
    
# Verticy coordinates  
  vert <- data.frame(x = c(0, -1/sqrt(2), 1/sqrt(2)),
                       y = c(sqrt(2)/sqrt(3), -sqrt(2)/(2*(sqrt(3))), -sqrt(2)/(2*(sqrt(3)))))
  
# Plot
  arg$x <- tridata$x
  arg$y <- tridata$y
  arg$xlab = ' '
  arg$ylab = ' '
  arg$bty = 'n'
  arg$axes = FALSE
  
  do.call(plot, arg)
  
# Add lines 
  segments(vert$x[1], vert$y[1], vert$x[2], vert$y[2], lwd = out.lwd, lty = out.lty, col = out.lcol)
  segments(vert$x[1], vert$y[1], vert$x[3], vert$y[3], lwd = out.lwd, lty = out.lty, col = out.lcol)
  segments(vert$x[2], vert$y[2], vert$x[3], vert$y[3], lwd = out.lwd, lty = out.lty, col = out.lcol)
  
# Origin
  if(isTRUE(achro)){
    points(x = 0, y = 0, pch = 15, col = achrocol, cex = achrosize)
  }
  
# Add text (coloured points better as in tcsplot?)
  if(isTRUE(labels)){
    text('S', x = -0.76, y = -0.39, xpd = TRUE, cex = cex.labels)
    text('M', x = 0, y = 0.88, xpd = TRUE, cex = cex.labels)
    text('L', x = 0.76, y = -0.39, xpd = TRUE, cex = cex.labels)
  }
  
}
```

#### Example 1: Dichromatism.

Reflectance data from four body regions of male and female *Ctenophorus ornatus* (Whiting et al. 2015, Biol J Linn Soc). Labium, throat, tongue, and mouth-roof.

**Q:** Which body regions are sexually dichromatic?

Calculate deltaS according to conspecific (tetrachromatic) visual system

``` r
specs <- list(lab = as.rspec(read.csv('data/dichromatism/lab.csv'), interp = FALSE),
              throat = as.rspec(read.csv('data/dichromatism/throat.csv'), interp = FALSE),
              roof = as.rspec(read.csv('data/dichromatism/roof.csv'), interp = FALSE),
              tongue = as.rspec(read.csv('data/dichromatism/tongue.csv'), interp = FALSE))
```

    ## wavelengths found in column 1 
    ## wavelengths found in column 1 
    ## wavelengths found in column 1 
    ## wavelengths found in column 1

``` r
# Ctenophorus ornatus
liz_vis <- sensmodel(c(360, 440, 493, 571)) 
names(liz_vis) <- c('wl', 'u', 's', 'm', 'l')

models <- lapply(specs, function(x) vismodel(x, visual = liz_vis, relative = FALSE, 
                                             qcatch = "fi", scale = 10000))  # deltaS
models_rel <- lapply(specs, function(x) vismodel(x, visual = liz_vis, relative = TRUE, 
                                                 qcatch = "fi", scale = 10000))  # tcs 

deltaS <- lapply(models, function(x) coldist(x, achro = FALSE, n1 = 1, n2 = 1, 
                                             n3 = 3.5, n4 = 6, v = 0.10))

# To add group labels (because I'm bad at R and I feel bad)
liz_lab <- function(x){
  x$comparison[grepl('F', x$patch1) & grepl('F', x$patch2)] <- 'intra.F'
  x$comparison[grepl('M', x$patch1) & grepl('M', x$patch2)] <- 'intra.M'
  x$comparison[grepl('M', x$patch1) & grepl('F', x$patch2)] <- 'inter'
  x$comparison[grepl('F', x$patch1) & grepl('M', x$patch2)] <- 'inter'
  x
}

# ew
deltaS$lab <- liz_lab(deltaS$lab) 
deltaS$throat <- liz_lab(deltaS$throat)
deltaS$roof <- liz_lab(deltaS$roof)
deltaS$tongue <- liz_lab(deltaS$tongue)
```

Plot 'em

``` r
par(pty="s", mfrow = c(2, 2))

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$lab[grepl("M", rownames(models_rel$lab)), ])
                                       [, c('x','y','z')]), pch=19, box=F, main = 'labium')
sp3d$points3d(suppressWarnings(tcs(models_rel$lab[grepl("F", rownames(models_rel$lab)), ])
                               [, c('x','y','z')]), col='red',pch=19)

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$throat[grepl("M", rownames(models_rel$throat)), ])
                                       [, c('x','y','z')]), pch=19, box=F, main = 'throat')
sp3d$points3d(suppressWarnings(tcs(models_rel$throat[grepl("F", rownames(models_rel$throat)), ])
                               [, c('x','y','z')]), col='red',pch=19)

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$roof[grepl("M", rownames(models_rel$roof)), ])
                                       [, c('x','y','z')]), pch=19, box=F, main = 'roof')
sp3d$points3d(suppressWarnings(tcs(models_rel$roof[grepl("F", rownames(models_rel$roof)), ])
                               [, c('x','y','z')]), col='red',pch=19)

sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel$tongue[grepl("M", rownames(models_rel$tongue)), ])
                                       [, c('x','y','z')]), pch=19, box=F, main = 'tongue')
sp3d$points3d(suppressWarnings(tcs(models_rel$tongue[grepl("F", rownames(models_rel$tongue)), ])
                               [, c('x','y','z')]), col='red',pch=19)
```

![](output/figures/examples/examples_figliz_tcs-1.png)

``` r
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

![](output/figures/examples/examples_figliz_deltaplot-1.png)

#### So what's sexually dichromatic?

**Ye olde ways**

Use mean inter-sex distance (eg. Bridge et al. 2008):

``` r
mean(deltaS$lab[deltaS$lab$comparison == 'inter', ][['dS']])
```

    ## [1] 5.059284

``` r
mean(deltaS$throat[deltaS$throat$comparison == 'inter', ][['dS']])
```

    ## [1] 6.083588

``` r
mean(deltaS$roof[deltaS$roof$comparison == 'inter', ][['dS']])
```

    ## [1] 2.514028

``` r
mean(deltaS$tongue[deltaS$tongue$comparison == 'inter', ][['dS']])
```

    ## [1] 3.025858

Or the maximums?? (eg. Burns & Schultz 2012):

``` r
max(deltaS$lab[deltaS$lab$comparison == 'inter', ][['dS']])
```

    ## [1] 12.71407

``` r
max(deltaS$throat[deltaS$throat$comparison == 'inter', ][['dS']])
```

    ## [1] 23.90123

``` r
max(deltaS$roof[deltaS$roof$comparison == 'inter', ][['dS']])
```

    ## [1] 13.76762

``` r
max(deltaS$tongue[deltaS$tongue$comparison == 'inter', ][['dS']])
```

    ## [1] 13.93182

Or run some t tests (eg. Igic et al. 2010):

``` r
t.test(deltaS$lab[deltaS$lab$comparison == 'inter', ][['dS']], mu = 1)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  deltaS$lab[deltaS$lab$comparison == "inter", ][["dS"]]
    ## t = 51.623, df = 863, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 1
    ## 95 percent confidence interval:
    ##  4.904951 5.213618
    ## sample estimates:
    ## mean of x 
    ##  5.059284

``` r
t.test(deltaS$throat[deltaS$throat$comparison == 'inter', ][['dS']], mu = 1)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  deltaS$throat[deltaS$throat$comparison == "inter", ][["dS"]]
    ## t = 42.348, df = 890, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 1
    ## 95 percent confidence interval:
    ##  5.847989 6.319187
    ## sample estimates:
    ## mean of x 
    ##  6.083588

``` r
t.test(deltaS$roof[deltaS$roof$comparison == 'inter', ][['dS']], mu = 1)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  deltaS$roof[deltaS$roof$comparison == "inter", ][["dS"]]
    ## t = 17.577, df = 809, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 1
    ## 95 percent confidence interval:
    ##  2.344950 2.683107
    ## sample estimates:
    ## mean of x 
    ##  2.514028

``` r
t.test(deltaS$tongue[deltaS$tongue$comparison == 'inter', ][['dS']], mu = 1)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  deltaS$tongue[deltaS$tongue$comparison == "inter", ][["dS"]]
    ## t = 24.356, df = 890, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 1
    ## 95 percent confidence interval:
    ##  2.862611 3.189105
    ## sample estimates:
    ## mean of x 
    ##  3.025858

...etc

**Conclusion**: Dimorphism everywhere!

#### New method

**Step 1:** Run permuational ANOVA (PERMANOVA) on lizard bits to ask if group A is different than group B

Labium

``` r
adoniscoldist(deltaS$lab)
```

    ## 
    ## Call:
    ## adonis(formula = dmat ~ grouping) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## grouping   1    150.12 150.119  14.117 0.20134  0.001 ***
    ## Residuals 56    595.50  10.634         0.79866           
    ## Total     57    745.62                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Throat

``` r
adoniscoldist(deltaS$throat)
```

    ## 
    ## Call:
    ## adonis(formula = dmat ~ grouping) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## grouping   1    202.58 202.583  14.978 0.20809  0.001 ***
    ## Residuals 57    770.97  13.526         0.79191           
    ## Total     58    973.55                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Mouth-roof

``` r
adoniscoldist(deltaS$roof)
```

    ## 
    ## Call:
    ## adonis(formula = dmat ~ grouping) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)
    ## grouping   1      3.22  3.2242 0.49025 0.009  0.502
    ## Residuals 54    355.14  6.5766         0.991       
    ## Total     55    358.36                 1.000

Tongue

``` r
adoniscoldist(deltaS$tongue)
```

    ## 
    ## Call:
    ## adonis(formula = dmat ~ grouping) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## grouping   1     12.17 12.1726  1.6766 0.02857  0.198
    ## Residuals 57    413.82  7.2601         0.97143       
    ## Total     58    426.00                 1.00000

**Conclusion**: labium = distinct, throat = distinct, mouth = nope, tongue = nope.

Much more sensible sensible based on the earlier plots & what the lizards actually look like (in the Whiting et al. paper). Neat.

**Step 2:** Effect sizes for distinct bits (i.e. do they crack the threshold?). Intercept-only mixed-models?

``` r
deltaS.lab.inter <- deltaS$lab[deltaS$lab$comparison == 'inter', ]

grep('M', deltaS.lab.inter$patch1)
```

    ## integer(0)

``` r
grep('F', deltaS.lab.inter$patch2)
```

    ## integer(0)

``` r
summary(lmer(dS~1+(1|patch1)+(1|patch2), data=deltaS.lab.inter))
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: dS ~ 1 + (1 | patch1) + (1 | patch2)
    ##    Data: deltaS.lab.inter
    ## 
    ## REML criterion at convergence: 3139.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.9358 -0.5955  0.0912  0.6083  2.4763 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  patch2   (Intercept) 0.7348   0.8572  
    ##  patch1   (Intercept) 2.9433   1.7156  
    ##  Residual             1.7919   1.3386  
    ## Number of obs: 864, groups:  patch2, 32; patch1, 27
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   5.0593     0.3661   13.82

**95% CI for labium = 4.69, 5.43**

``` r
deltaS.throat.inter <- deltaS$throat[deltaS$throat$comparison == 'inter', ]

grep('M', deltaS.throat.inter$patch1)
```

    ## integer(0)

``` r
grep('F', deltaS.throat.inter$patch2)
```

    ## integer(0)

``` r
summary(lmer(dS~1+(1|patch1)+(1|patch2), data=deltaS.throat.inter))
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: dS ~ 1 + (1 | patch1) + (1 | patch2)
    ##    Data: deltaS.throat.inter
    ## 
    ## REML criterion at convergence: 3393
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.6384 -0.5730  0.0502  0.6138  3.5456 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  patch2   (Intercept) 7.797    2.792   
    ##  patch1   (Intercept) 3.419    1.849   
    ##  Residual             1.974    1.405   
    ## Number of obs: 891, groups:  patch2, 33; patch1, 27
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   6.0836     0.6043   10.07

**95% CI for throat = 5.48, 6.69**

``` r
rm(deltaS, deltaS.lab.inter, deltaS.throat.inter, models, models_rel, specs, liz_vis, liz_lab)
```

#### Example 2: Mimicry.

Reflectance data from colour-polymorphic female spiders *Gasteracantha fornicata*, and sympatic flowers from Qld, Australia. (W = white morph, Y = yellow morph, F = flowers)

So three groups, with **two Q's:**

**(1)** Are spiders actually polymorphic (to prey)?

**(2)** Do spiders (of each morph) resemble sympatric flowers?

``` r
# specs <- list(spider_w = as.rspec(read.csv('data/mimicry/spiders_white.csv'), interp = FALSE),
#               spider_y = as.rspec(read.csv('data/mimicry/spiders_yellow.csv'), interp = FALSE),
#               flower = as.rspec(read.csv('data/mimicry/flowers.csv'), interp = FALSE))

specs <- as.rspec(read.csv('data/mimicry/flowers_spiders.csv'), interp = FALSE)
```

    ## wavelengths found in column 1

Calculate deltaS according to bee (trichromat) visual system

``` r
# Honeybee
bee_vis <- sensmodel(c(350, 440, 540)) 
names(bee_vis) <- c('wl', 's', 'm', 'l')

models <- vismodel(specs, visual = bee_vis, relative = FALSE, 
                                             qcatch = "fi", scale = 10000)  # deltaS
models_rel <- vismodel(specs, visual = bee_vis, relative = TRUE, 
                                                 qcatch = "fi", scale = 10000)  # max triangle
models_tri <- trispace(models_rel)

deltaS <- coldist(models, achro = FALSE)

# Contrast labels
deltaS$comparison[grepl('W_', deltaS$patch1) & grepl('W_', deltaS$patch2)] <- 'intra.W'
deltaS$comparison[grepl('Y_', deltaS$patch1) & grepl('Y_', deltaS$patch2)] <- 'intra.Y'
deltaS$comparison[grepl('F_', deltaS$patch1) & grepl('F_', deltaS$patch2)] <- 'intra.F'
deltaS$comparison[grepl('W_', deltaS$patch1) & grepl('Y_', deltaS$patch2)] <- 'inter.WY'
deltaS$comparison[grepl('Y_', deltaS$patch1) & grepl('W_', deltaS$patch2)] <- 'inter.WY'
deltaS$comparison[grepl('W_', deltaS$patch1) & grepl('F_', deltaS$patch2)] <- 'inter.WF'
deltaS$comparison[grepl('F_', deltaS$patch1) & grepl('W_', deltaS$patch2)] <- 'inter.WF'
deltaS$comparison[grepl('F_', deltaS$patch1) & grepl('Y_', deltaS$patch2)] <- 'inter.YF'
deltaS$comparison[grepl('Y_', deltaS$patch1) & grepl('F_', deltaS$patch2)] <- 'inter.YF'
```

``` r
triplot(models_tri[grepl("W_", rownames(models_tri)), ], col = 'darkgrey')
points(subset(models_tri[grepl("Y_", rownames(models_tri)), ], select = c('x', 'y')), pch = 19, col = 'yellow')
points(subset(models_tri[grepl("F_", rownames(models_tri)), ], select = c('x', 'y')), pch = 19, col = 'forestgreen')
```

![](output/figures/examples/examples_figmimic_triplot-1.png)

Well that's not terribly helpful. Should actually run it in the hexagon w/ euclidean distances since it's a bee model.

``` r
ggplot(deltaS, aes(x=dS, fill=comparison)) + geom_histogram(bins=50) + 
        facet_grid(comparison~., scales='free_y') + geom_vline(xintercept=1) +
        ggtitle('labial') + theme(legend.position="none")
```

![](output/figures/examples/examples_figmimic_deltaplot-1.png)

Mmm, maybe not the most exciting example.

...tbc

``` r
rm(deltaS, models, models_rel, models_tri, specs, bee_vis)
```

#### Example 3: Crypsis.

Reflectance data from various body regions (H = head, L = left arm, R = right arm, P = prothorax, W = wing, A = abdomen) of 27 female mantids *Pseudomantis albofimbriata* and 50 background samples (*Lomandra longifolia*, which they pretty much exclusively hang on).

So six groups, one **Q:** Are mantids cryptic? i.e. are all body regions chromaticically indistinguishable from their background?

Calculate deltaS according to blue tits

``` r
specs <- as.rspec(read.csv('data/crypsis/mantids_bkgs.csv'), lim = c(300, 700))
```

    ## wavelengths found in column 1

``` r
models <- vismodel(specs, visual = 'bluetit', relative = FALSE, qcatch = "fi", scale = 10000)  # deltaS
models_rel <- vismodel(specs, visual = 'bluetit', relative = TRUE, qcatch = "fi", scale = 10000, achro = FALSE)  # tcs

deltaS <- coldist(models, achro = FALSE)
```

Plot 'em

``` r
sp3d <- scatterplot3d(suppressWarnings(tcs(models_rel[grepl("B_", rownames(models_rel)), ])
                                       [, c('x','y','z')]), xlim=c(-0.02,0.02), ylim=c(-0.01,0.028), zlim=c(-0.06,0.01), 
                                        pch=19, box=F, color = 'forestgreen')
sp3d$points3d(suppressWarnings(tcs(models_rel[grepl("H_", rownames(models_rel)), ])
                               [, c('x','y','z')]), col='darkgoldenrod',pch=19)
sp3d$points3d(suppressWarnings(tcs(models_rel[grepl("L_", rownames(models_rel)), ])
                               [, c('x','y','z')]), col='darkgoldenrod1',pch=19)
sp3d$points3d(suppressWarnings(tcs(models_rel[grepl("R_", rownames(models_rel)), ])
                               [, c('x','y','z')]), col='darkgoldenrod2',pch=19)
sp3d$points3d(suppressWarnings(tcs(models_rel[grepl("P_", rownames(models_rel)), ])
                               [, c('x','y','z')]), col='darkgoldenrod3',pch=19)
sp3d$points3d(suppressWarnings(tcs(models_rel[grepl("W_", rownames(models_rel)), ])
                               [, c('x','y','z')]), col='darkgoldenrod4',pch=19)
sp3d$points3d(suppressWarnings(tcs(models_rel[grepl("A_", rownames(models_rel)), ])
                               [, c('x','y','z')]), col='gold1',pch=19)
```

![](output/figures/examples/examples_figmantid_tcs-1.png)

...tbc

``` r
sessionInfo()
```

    ## R version 3.3.1 (2016-06-21)
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
    ##  [1] lme4_1.1-12          Matrix_1.2-6         vegan_2.4-0         
    ##  [4] lattice_0.20-33      permute_0.9-0        gridExtra_2.2.1     
    ##  [7] ggplot2_2.1.0        scatterplot3d_0.3-37 pavo_0.5-5          
    ## [10] rgl_0.95.1441       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.5      nloptr_1.0.4     formatR_1.4      plyr_1.8.4      
    ##  [5] tools_3.3.1      magic_1.5-6      digest_0.6.9     evaluate_0.9    
    ##  [9] gtable_0.2.0     nlme_3.1-128     mgcv_1.8-12      mapproj_1.2-4   
    ## [13] yaml_2.1.13      parallel_3.3.1   stringr_1.0.0    knitr_1.13      
    ## [17] cluster_2.0.4    maps_3.1.0       rcdd_1.1-10      grid_3.3.1      
    ## [21] rmarkdown_0.9.6  minqa_1.2.4      reshape2_1.4.1   magrittr_1.5    
    ## [25] scales_0.4.0     htmltools_0.3.5  MASS_7.3-45      splines_3.3.1   
    ## [29] colorspace_1.2-6 labeling_0.3     stringi_1.1.1    geometry_0.3-6  
    ## [33] munsell_0.4.3
