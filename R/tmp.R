
Calculate statistics of interest:
  
  ```{r, dependson=c('powercoldistcache', 'poweradoniscache', 'powerpyke', 'powercoldistcache2', 'poweradoniscache2', 'powerpyke2')}
# Pyke MANOVA: P-values, which significant
manovaPval.n10 <- unlist(lapply(pykemanova.n10, function(x) x$stats[1,'Pr(>F)']))
manovaPval.n20 <- unlist(lapply(pykemanova.n20, function(x) x$stats[1,'Pr(>F)']))
manovaPval.n100 <- unlist(lapply(pykemanova.n100, function(x) x$stats[1,'Pr(>F)']))
manovaP.n10 <- manovaPval.n10 < 0.05
manovaP.n20 <- manovaPval.n20 < 0.05
manovaP.n100 <- manovaPval.n100 < 0.05


# Adonis: P-values, which significant, R-squared
adonisPval.n10 <- unlist(lapply(adonissim.n10, function(x) x$aov.tab$'Pr(>F)'[1]))
adonisPval.n20 <- unlist(lapply(adonissim.n20, function(x) x$aov.tab$'Pr(>F)'[1]))
adonisPval.n100 <- unlist(lapply(adonissim.n100, function(x) x$aov.tab$'Pr(>F)'[1]))
adonisP.n10 <- adonisPval.n10 < 0.05
adonisP.n20 <- adonisPval.n20 < 0.05
adonisP.n100 <- adonisPval.n100 < 0.05
adonisR2.n10 <- unlist(lapply(adonissim.n10, function(x) x$aov.tab$'R2'[1])) * 100
adonisR2.n20 <- unlist(lapply(adonissim.n20, function(x) x$aov.tab$'R2'[1])) * 100
adonisR2.n100 <- unlist(lapply(adonissim.n100, function(x) x$aov.tab$'R2'[1])) * 100

# Disagreement in significance between Pyke and Adonis
disagreement.n10 <- !manovaP.n10 == adonisP.n10
disagreement.n20 <- !manovaP.n20 == adonisP.n20
disagreement.n100 <- !manovaP.n100 == adonisP.n100

# Proportion of tests significant by effect size, according to each method
pvsimsadonis.n10 <- tapply(adonisP.n10, effsims, mean)
pvsimsadonis.n20 <- tapply(adonisP.n20, effsims, mean)
pvsimsadonis.n100 <- tapply(adonisP.n100, effsims, mean)
pvsimsmanova.n10 <- tapply(manovaP.n10, effsims, mean)
pvsimsmanova.n20 <- tapply(manovaP.n20, effsims, mean)
pvsimsmanova.n100 <- tapply(manovaP.n100, effsims, mean)

```

##Visualizing Results

```{r, echo=FALSE, dependson=c('powercoldistcache', 'poweradoniscache', 'powerpyke')}
powerpal <- brewer.pal(9,'PRGn')

plot(pvsimsadonis.n10~unique(effsims), pch=20, cex=2, type='b', ylab='Proportion of tests significant', xlab="Effect size (Mahalanobis' Distance)", ylim=c(0,1), col=powerpal[6])

points(pvsimsadonis.n20~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[7])
points(pvsimsadonis~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[8])
points(pvsimsadonis.n100~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[9])

points(pvsimsmanova.n10~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[4])
points(pvsimsmanova.n20~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[3])
points(pvsimsmanova~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[2])
points(pvsimsmanova.n100~ unique(effsims), pch=20, cex=2, type='b', col=powerpal[1])

legend('topleft', legend=c("adonis, N=10",
                           "adonis, N=20",
                           "adonis, N=50",
                           "adonis, N=100",
                           "Pyke & MANOVA, N=10",
                           "Pyke & MANOVA, N=20",
                           "Pyke & MANOVA, N=50",
                           "Pyke & MANOVA, N=100"), 
       pch=20, pt.cex=2, lty=1, bty='n', seg.len=3,
       col=powerpal[c(6,7,8,9,4,3,2,1)])
``` 
