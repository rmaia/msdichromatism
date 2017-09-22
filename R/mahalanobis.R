mahalanobis <- function (data, groups, mve = TRUE) {
	
	groups <- factor(groups)

    meanA <- colMeans(data[groups %in% levels(groups)[1],])
    meanB <- colMeans(data[groups %in% levels(groups)[2],])
    
    
    if (mve) {
    	nA <- dim(data[groups %in% levels(groups)[1],])[1]
    	nB <- dim(data[groups %in% levels(groups)[2],])[1]
    	
        covarA <- cov.mve(data[groups %in% levels(groups)[1],])$cov
        covarB <- cov.mve(data[groups %in% levels(groups)[2],])$cov
        covarpool <- ((nA-1)*covarA + (nB-1)*covarB)/(nA+nB-1)
    }
    else {
    	nA <- dim(data[groups %in% levels(groups)[1],])[1]
    	nB <- dim(data[groups %in% levels(groups)[2],])[1]
    	
        covarA <- var(data[groups %in% levels(groups)[1],])
        covarB <- var(data[groups %in% levels(groups)[2],])
        covarpool <- ((nA-1)*covarA + (nB-1)*covarB)/(nA+nB-1)
    }

    t(meanA - meanB) %*% solve(covarpool) %*% (meanA - meanB)
}