# EXAMPLE AT THE END OF THE SCRIPT
# dat = vismodel object (quantum catch)
# groups = vector defining the groups to be compared
# nboot = number of bootstrap replicates
# alpha = confidence interval
# ... = attributes passed to coldist

bootcentroidDS <- function(dat, groups, nboot=1000, alpha=0.95, ...){

sortinggroups <- order(groups)
dat <- dat[sortinggroups, ]
groups <- groups[sortinggroups]
  
  
samplesizes <- table(groups)

# calculate empirical deltaS
empgroupmeans <- aggregate(dat, list(groups), mean, simplify=TRUE)
row.names(empgroupmeans) <- empgroupmeans[,1]
empgroupmeans <- empgroupmeans[, -1]

empcd <- coldist(empgroupmeans, ...)
empdS <- setNames(empcd$dS, paste(empcd$patch1,empcd$patch2, sep='-'))


# separate data by group
bygroup <- lapply(unique(groups), function(x) dat[groups==x, ])
# split(dat, groups) also works but is about twice as slow
names(bygroup) <- unique(groups)

# create vectors of indices to sample
its <- lapply(samplesizes, function(x) sample(1:x, x*nboot, replace=TRUE))

# sample
# returns a list with entries = number of groups
# and rows equal to the sample size for that group times the number of bootstrap replicates
bootsamples <- lapply(1:length(bygroup), function(x) bygroup[[x]][its[[x]], ] )
#bootsamples <- lapply(1:length(bygroup), function(x) bygroup[[names(its[x])]][its[[x]], ] )
# next we need to split by bootstrap replicate
# preserving the same sample size as that original group had

# first dimension of the list is the original group
# second dimension is the bootstrap replicate
bootindex <- lapply(samplesizes, function(x) as.character(rep(1:nboot, each=x)))

bootbygroup <- lapply(1:length(bygroup), function(x){
  lapply(unique(bootindex[[x]]), function(z) bootsamples[[x]][bootindex[[x]]==z, ])  
}
)

# should work but isn't:
#bootbygroup <- lapply(1:length(bygroup), function(x) 
#  split(bootsamples[[x]], rep(1:nboot, each=samplesizes[x] )
#       ))




# now we need to take the column means for all of these
groupcolmeans <- lapply(bootbygroup, function(x) do.call(rbind, lapply(x, colMeans) ) )

# ...and combine them by bootstrap replicate
bootgrouped <- lapply(1:nboot, function(x) 
  do.call(rbind, lapply(groupcolmeans, '[', x, )))

# ...name the rows by group
bootgrouped <- lapply(bootgrouped, function(x){row.names(x) <- unique(groups); x})

# apply coldist to the replicated bootstraps
bootcd <- lapply(bootgrouped, coldist, ...)

# get deltaS and name by group difference
bootdS <- do.call(rbind,
                  lapply(bootcd, function(x) 
                  setNames(x$dS, paste(x$patch1, x$patch2, sep='-')))
                  )

# ...subtract them from the empirical
# bootdS <- bootdS - empdS

# ... sort and find quantiles
#bootdS <- sort(bootdS)
#quantileindices <- round(length(bootdS)*((1+c(-alpha, alpha))/2))
#dsCI <- empdS + bootdS[quantileindices]

# which gives the same as just
quantileindices <- round(nboot*((1+c(-alpha, alpha))/2))
#bootdS <- sort(bootdS)
bootdS <- apply(bootdS, 2, sort)
dsCI <- bootdS[quantileindices, , drop=FALSE]
rownames(dsCI) <- c('CI.lwr','CI.upr')
# make sure names match with empirical (they always should but just in case)
dsCI <- dsCI[, names(empdS), drop=FALSE]

measured.dS <- empdS

t(rbind(measured.dS, dsCI))

}


# EXAMPLE
# (commented so function can be sourced)

# for example where groups aren't different, uncomment this
# test <- matrix(rnorm(400, 10, 3), ncol=4); colnames(test) <- c('u','s','m','l')

# for example where groups are different, uncomment this
#test <- rbind( matrix(rnorm(200, 10), ncol=4), matrix(c(rnorm(100, 5), rnorm(100, 10)), ncol=4)); colnames(test) <- c('u','s','m','l')

# for example with two groups
#grps <- rep(c('ga','gb'), each=50)

# for example with four groups
#grps <- rep(c('ga','gb', 'gc', 'gd'), each=25)

 #bootcentroidDS(test, grps)
# bootcentroidDS(test, grps, v=0.4)