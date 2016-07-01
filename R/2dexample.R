# illustrating the problem in two dimensions

png('figures/2dexample.png')
set.seed(9756)

movePoints <- function(x0y0, xy, d){
  total.dist <- apply(cbind(x0y0, xy), 1,
                      function(x) stats::dist(rbind(x[1:2], x[3:4])))
  p <- d / total.dist
  p <- 1 - p
  x0y0[,1] <- xy[,1] + p*(x0y0[,1] - xy[,1])
  x0y0[,2] <- xy[,2] + p*(x0y0[,2] - xy[,2])
  return(x0y0)
}

gA <- MASS::mvrnorm(10, mu=c(10,10), Sigma=matrix(c(1,0,0,1), nrow=2))
gB <- MASS::mvrnorm(10, mu=c(10,10), Sigma=matrix(c(1,0,0,1), nrow=2))

gB2 <- movePoints(gB, matrix(rep(gA[1,], each=10), ncol=2), 0.1)

par(pty='s')
par(mar=c(1,1,1,1))

plot(1,1, type='n', ylim=c(7.5,12), xlim=c(7.5, 11),
     yaxt='n', xaxt='n', ylab='', xlab='')

arrows(gA[1,1], gA[1,2],
       gB2[,1], gB2[,2],
       length=0.1, angle=30, lwd=1,
       col=grey(0.2, 0.8))

points(rbind(gA,gB), col=rep(c('#4575b4','#d73027'), each=10), 
       pch=19, cex=1.5)

points(matrix(colMeans(gA), ncol=2), pch=3, lwd=3, col='#4575b4', cex=2)
points(matrix(colMeans(gB), ncol=2), pch=3, lwd=3, col='#d73027', cex=2)

arrows( x0=colMeans(gA)[1], y0=colMeans(gA)[2], 
        x1=colMeans(gB)[1], y1=colMeans(gB)[2],
  length=0.1, angle=20, lwd=2, code=3, col=grey(0,0.8)
)

dev.off()