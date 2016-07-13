
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
  
  # Descriptive attributes (largely preserved from vismodel)
  # attr(res, 'clrsp') <- 'trispace'
  # attr(res, 'conenumb') <- 3
  # attr(res, 'qcatch') <- attr(vismodeldata, 'qcatch')
  # attr(res,'visualsystem.chromatic') <- attr(vismodeldata,'visualsystem.chromatic')
  # attr(res,'visualsystem.achromatic') <- attr(vismodeldata,'visualsystem.achromatic')
  # attr(res,'illuminant') <- attr(vismodeldata,'illuminant')
  # attr(res,'background') <- attr(vismodeldata,'background')
  # attr(res,'relative') <- attr(vismodeldata,'relative')
  # attr(res, 'vonkries') <- attr(vismodeldata, 'vonkries')
  # 
  # # Data attributes
  # attr(res, 'data.visualsystem.chromatic') <- attr(vismodeldata, 'data.visualsystem.chromatic')
  # attr(res, 'data.visualsystem.achromatic') <- attr(vismodeldata, 'data.visualsystem.achromatic')
  # attr(res, 'data.background') <- attr(vismodeldata, 'data.background')
  
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
