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
require(vegan)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.4-0

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
set.seed(1982)

# we'll just consider usml are uncorrelated for this example

# we'll generate data from a lognormal distribution to avoid negative values
# now with a variance (on the log scale) of 0.02

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

Again, although USML were simulated uncorrelated, XYZ are correlated.

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

``` r
deltaS.nodiff <- deltaS
alldat.nodiff <- alldat
```

#### Note that in both cases there is no dichromatism, but regular analyses would consider the second case dichromatic - because mean deltaS between males and females is over the threshold of 1 JND.

### Example 3: low intra-group variability, high inter-group distance

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
  m = rlnorm(30, meanlog=log(0.7), sdlog=sqrt(0.002)),
  l = rlnorm(30, meanlog=log(0.3), sdlog=sqrt(0.002))
)


colnames(groupA) <- colnames(groupB) <- c('u','s','m', 'l')
attr(groupA, 'relative') <- attr(groupB, 'relative') <- FALSE

par(pty="s")
sp3d <- scatterplot3d(suppressWarnings(tcs(groupA)[, c('x','y','z')]), pch=19,
                      xlim=c(0.1,0.38), ylim=c(-0.1,0.32), zlim=c(-0.2,-0.16), box=F)
sp3d$points3d(suppressWarnings(tcs(groupB)[, c('x','y','z')]), col='red',pch=19)
```

![](Simulations_files/figure-markdown_github/unnamed-chunk-6-1.png)

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

![](Simulations_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
deltaS.diff <- deltaS
alldat.diff <- alldat
```

How to test this?
-----------------

I think something from the likes of the adonis function in vegan, but I don't know how to use a custom distance matrix with it yet.

Really rough possibility using mixed models and crossed random effects:

#### when there's no difference

``` r
require(lme4)
```

    ## Loading required package: lme4

    ## Loading required package: Matrix

``` r
# no gain in explanatory power
anova(
  lmer(dS~(1|patch1)+(1|patch2), data=deltaS.nodiff),
  lmer(dS~comparison+(1|patch1)+(1|patch2), data=deltaS.nodiff),
  model.names = c('null','alternative')
)
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: deltaS.nodiff
    ## Models:
    ## null: dS ~ (1 | patch1) + (1 | patch2)
    ## alternative: dS ~ comparison + (1 | patch1) + (1 | patch2)
    ##             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## null         4 4182.4 4204.3 -2087.2   4174.4                         
    ## alternative  6 4183.2 4216.1 -2085.6   4171.2 3.1536      2     0.2066

``` r
summary(lmer(dS~comparison-1+(1|patch1)+(1|patch2), data=deltaS.nodiff))
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: dS ~ comparison - 1 + (1 | patch1) + (1 | patch2)
    ##    Data: deltaS.nodiff
    ## 
    ## REML criterion at convergence: 4178.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4491 -0.5360  0.0786  0.7126  2.3055 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  patch1   (Intercept) 0.2484   0.4984  
    ##  patch2   (Intercept) 0.1875   0.4330  
    ##  Residual             0.5289   0.7272  
    ## Number of obs: 1770, groups:  patch1, 59; patch2, 59
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error t value
    ## comparisoninter     2.1695     0.1229   17.65
    ## comparisonintra.A   2.3254     0.1298   17.91
    ## comparisonintra.B   2.0016     0.1304   15.35
    ## 
    ## Correlation of Fixed Effects:
    ##             cmprsn cmpr.A
    ## cmprsnntr.A 0.519        
    ## cmprsnntr.B 0.390  0.000

According to Eaton 2005:

> Because average reflectance curves were used in the color discrimination model, between-sex differences identified by the model might not be biologically functional if variance in coloration within sexes is so broad as not to be a reliable visual indicator of sex. Hence, I assessed intraspecific variation in coloration between sexes using logistic regression (PROC GENMOD, SAS V.8.0, SAS Institute, Cary, NC), with sex (1 = male, 0 = female) as the response variable and Qi (i.e., receptor quantum catches, Eq. 1 ) as predictor variables. I modeled the probability of an individual being male given a value for Qi , for each of the four receptor quantum catches for each feather patch within each species. If the model regression coefficient estimate was zero, then that quantum catch had no effect on sex (i.e., it cannot predict sex). A positive regression coefficient indicated an increased probability of an individual being male with larger values of Qi , whereas a negative regression coefficient indicated a higher probability of being female with larger values of Qi . I used likelihood ratio confidence intervals for estimating whether strong correlations existed between the response variable (sex) and the predictor variables (Q1–Q4) for each feather patch. Given the small sample sizes for each feather patch comparison (n = 10), I report 85% upper and lower confidence intervals around the regression coefficient estimates (29) (see Table 1).

``` r
alldat.nodiff$group <- gsub('[0-9]', '',rownames(alldat.nodiff))

summary(glm(as.numeric(as.factor(group))~u, data=alldat.nodiff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ u, data = alldat.nodiff)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.77670  -0.46423  -0.04167   0.43751   0.71222  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   0.5317     0.5234   1.016   0.3140  
    ## u             9.5493     5.1236   1.864   0.0674 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.244007)
    ## 
    ##     Null deviance: 15.000  on 59  degrees of freedom
    ## Residual deviance: 14.152  on 58  degrees of freedom
    ## AIC: 89.605
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(as.numeric(as.factor(group))~s, data=alldat.nodiff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ s, data = alldat.nodiff)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.51060  -0.50101  -0.00176   0.50098   0.50819  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)   1.4703     0.4682   3.141  0.00265 **
    ## s             0.2862     4.4628   0.064  0.94909   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.2586024)
    ## 
    ##     Null deviance: 15.000  on 59  degrees of freedom
    ## Residual deviance: 14.999  on 58  degrees of freedom
    ## AIC: 93.091
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(as.numeric(as.factor(group))~m, data=alldat.nodiff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ m, data = alldat.nodiff)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.59657  -0.49701   0.01609   0.49458   0.54262  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.7835     0.4827   3.695 0.000489 ***
    ## m            -0.9396     1.5848  -0.593 0.555562    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.2570627)
    ## 
    ##     Null deviance: 15.00  on 59  degrees of freedom
    ## Residual deviance: 14.91  on 58  degrees of freedom
    ## AIC: 92.732
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(as.numeric(as.factor(group))~l, data=alldat.nodiff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ l, data = alldat.nodiff)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.62230  -0.47897  -0.02255   0.50572   0.63300  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   1.0188     0.4859   2.097   0.0404 *
    ## l             0.6649     0.6654   0.999   0.3218  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.2542437)
    ## 
    ##     Null deviance: 15.000  on 59  degrees of freedom
    ## Residual deviance: 14.746  on 58  degrees of freedom
    ## AIC: 92.071
    ## 
    ## Number of Fisher Scoring iterations: 2

Seems like this would work but does not consider USML will be correlated (which they always will) and probably subject to type I error overinflation.

``` r
dmat <- matrix(0, nrow=length(unique(deltaS.nodiff$patch1)), ncol=length(unique(deltaS.nodiff$patch1)))
rownames(dmat) <- colnames(dmat) <- as.character(unique(deltaS.nodiff$patch1))

for(i in rownames(dmat))
  for(j in colnames(dmat))
    if(length(deltaS.nodiff$dS[deltaS.nodiff$patch1 == i & deltaS.nodiff$patch2 == j]) != 0)
      dmat[i,j] <- dmat[j,i] <- deltaS.nodiff$dS[deltaS.nodiff$patch1 == i & deltaS.nodiff$patch2 == j]

grouping <- gsub('[0-9]','', rownames(dmat))

anosim(dmat, grouping)
```

    ## 
    ## Call:
    ## anosim(dat = dmat, grouping = grouping) 
    ## Dissimilarity: user supplied square matrix 
    ## 
    ## ANOSIM statistic R: 0.001287 
    ##       Significance: 0.394 
    ## 
    ## Permutation: free
    ## Number of permutations: 999

#### when there is a difference

``` r
# gain in explanatory power
anova(
  lmer(dS~(1|patch1)+(1|patch2), data=deltaS.diff),
  lmer(dS~comparison+(1|patch1)+(1|patch2), data=deltaS.diff),
  model.names = c('null','alternative')
)
```

    ## refitting model(s) with ML (instead of REML)

    ## Data: deltaS.diff
    ## Models:
    ## null: dS ~ (1 | patch1) + (1 | patch2)
    ## alternative: dS ~ comparison + (1 | patch1) + (1 | patch2)
    ##             Df     AIC     BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
    ## null         4 1175.59 1197.50 -583.79  1167.59                         
    ## alternative  6  500.66  533.53 -244.33   488.66 678.93      2  < 2.2e-16
    ##                
    ## null           
    ## alternative ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(lmer(dS~comparison-1+(1|patch1)+(1|patch2), data=deltaS.diff))
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: dS ~ comparison - 1 + (1 | patch1) + (1 | patch2)
    ##    Data: deltaS.diff
    ## 
    ## REML criterion at convergence: 499.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.3848 -0.5138  0.0016  0.4987  4.6611 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  patch1   (Intercept) 0.07800  0.2793  
    ##  patch2   (Intercept) 0.05311  0.2305  
    ##  Residual             0.06246  0.2499  
    ## Number of obs: 1770, groups:  patch1, 59; patch2, 59
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error t value
    ## comparisoninter    9.92017    0.06663  148.88
    ## comparisonintra.A  0.66347    0.06865    9.66
    ## comparisonintra.B  0.54482    0.06894    7.90
    ## 
    ## Correlation of Fixed Effects:
    ##             cmprsn cmpr.A
    ## cmprsnntr.A 0.568        
    ## cmprsnntr.B 0.385  0.000

According to Eaton 2005:

> Because average reflectance curves were used in the color discrimination model, between-sex differences identified by the model might not be biologically functional if variance in coloration within sexes is so broad as not to be a reliable visual indicator of sex. Hence, I assessed intraspecific variation in coloration between sexes using logistic regression (PROC GENMOD, SAS V.8.0, SAS Institute, Cary, NC), with sex (1 = male, 0 = female) as the response variable and Qi (i.e., receptor quantum catches, Eq. 1 ) as predictor variables. I modeled the probability of an individual being male given a value for Qi , for each of the four receptor quantum catches for each feather patch within each species. If the model regression coefficient estimate was zero, then that quantum catch had no effect on sex (i.e., it cannot predict sex). A positive regression coefficient indicated an increased probability of an individual being male with larger values of Qi , whereas a negative regression coefficient indicated a higher probability of being female with larger values of Qi . I used likelihood ratio confidence intervals for estimating whether strong correlations existed between the response variable (sex) and the predictor variables (Q1–Q4) for each feather patch. Given the small sample sizes for each feather patch comparison (n = 10), I report 85% upper and lower confidence intervals around the regression coefficient estimates (29) (see Table 1).

``` r
alldat.diff$group <- gsub('[0-9]', '',rownames(alldat.diff))

summary(glm(as.numeric(as.factor(group))~u, data=alldat.diff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ u, data = alldat.diff)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.75406  -0.46923  -0.02235   0.43497   0.72179  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)   -1.507      1.650  -0.914   0.3647  
    ## u             29.989     16.440   1.824   0.0733 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.2445883)
    ## 
    ##     Null deviance: 15.000  on 59  degrees of freedom
    ## Residual deviance: 14.186  on 58  degrees of freedom
    ## AIC: 89.748
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(as.numeric(as.factor(group))~s, data=alldat.diff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ s, data = alldat.diff)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.5089  -0.5011  -0.0011   0.5007   0.5081  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   1.4145     1.4684   0.963    0.339
    ## s             0.8465    14.5257   0.058    0.954
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.2586055)
    ## 
    ##     Null deviance: 15.000  on 59  degrees of freedom
    ## Residual deviance: 14.999  on 58  degrees of freedom
    ## AIC: 93.091
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(as.numeric(as.factor(group))~m, data=alldat.diff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ m, data = alldat.diff)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -0.082289  -0.041081  -0.000665   0.033397   0.122495  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.25523    0.01698   15.03   <2e-16 ***
    ## m            2.49240    0.03156   78.97   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.002383398)
    ## 
    ##     Null deviance: 15.00000  on 59  degrees of freedom
    ## Residual deviance:  0.13824  on 58  degrees of freedom
    ## AIC: -188.12
    ## 
    ## Number of Fisher Scoring iterations: 2

``` r
summary(glm(as.numeric(as.factor(group))~l, data=alldat.diff))
```

    ## 
    ## Call:
    ## glm(formula = as.numeric(as.factor(group)) ~ l, data = alldat.diff)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -0.137118  -0.033269  -0.001287   0.032821   0.137298  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2.74956    0.01888  145.65   <2e-16 ***
    ## l           -2.48297    0.03485  -71.24   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.002922144)
    ## 
    ##     Null deviance: 15.00000  on 59  degrees of freedom
    ## Residual deviance:  0.16948  on 58  degrees of freedom
    ## AIC: -175.89
    ## 
    ## Number of Fisher Scoring iterations: 2

Seems like this would work but does not consider USML will be correlated (which they always will) and probably subject to type I error overinflation.

``` r
dmat <- matrix(0, nrow=length(unique(deltaS.diff$patch1)), ncol=length(unique(deltaS.diff$patch1)))
rownames(dmat) <- colnames(dmat) <- as.character(unique(deltaS.diff$patch1))

for(i in rownames(dmat))
  for(j in colnames(dmat))
    if(length(deltaS.diff$dS[deltaS.diff$patch1 == i & deltaS.diff$patch2 == j]) != 0)
      dmat[i,j] <- dmat[j,i] <- deltaS.diff$dS[deltaS.diff$patch1 == i & deltaS.diff$patch2 == j]

grouping <- gsub('[0-9]','', rownames(dmat))

anosim(dmat, grouping)
```

    ## 
    ## Call:
    ## anosim(dat = dmat, grouping = grouping) 
    ## Dissimilarity: user supplied square matrix 
    ## 
    ## ANOSIM statistic R:     1 
    ##       Significance: 0.001 
    ## 
    ## Permutation: free
    ## Number of permutations: 999

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
    ## [1] lme4_1.1-12          Matrix_1.2-6         vegan_2.4-0         
    ## [4] lattice_0.20-33      permute_0.9-0        ggplot2_2.1.0       
    ## [7] scatterplot3d_0.3-37 pavo_0.5-5           rgl_0.95.1441       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.5        nloptr_1.0.4       plyr_1.8.3        
    ##  [4] tools_3.3.0        magic_1.5-6        digest_0.6.9      
    ##  [7] evaluate_0.9       gtable_0.2.0       nlme_3.1-127      
    ## [10] mgcv_1.8-12        mapproj_1.2-4      yaml_2.1.13       
    ## [13] parallel_3.3.0     stringr_1.0.0      knitr_1.13        
    ## [16] cluster_2.0.4      maps_3.1.0         rcdd_1.1-10       
    ## [19] grid_3.3.0         rmarkdown_0.9.6.10 minqa_1.2.4       
    ## [22] reshape2_1.4.1     magrittr_1.5       scales_0.4.0      
    ## [25] htmltools_0.3.5    MASS_7.3-45        splines_3.3.0     
    ## [28] colorspace_1.2-6   labeling_0.3       stringi_1.0-1     
    ## [31] geometry_0.3-6     munsell_0.4.3
