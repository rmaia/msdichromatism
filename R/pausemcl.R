pausemcl <- function(X, FUN, splt=5, mcr=4){
  N <- length(X)
  resultvector <- vector('list', length=N)
  splitlist <- parallel::splitIndices(N, splt)
  
  for(i in seq_along(splitlist)){
    tmp <- splitlist[[i]]
    gc()
    resultvector[tmp] <- parallel::mclapply(X[tmp], FUN, mc.cores=mcr)
    #resultvector[tmp] <- lapply(X[tmp], FUN)
    gc()
    Sys.sleep(2)
  }
  resultvector
}
