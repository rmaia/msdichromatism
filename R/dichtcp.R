dichtcp <- function(dat, angle=70, scale.y=0.45){
  
  oldpar.pty <- par()$pty
  oldpar.mar <- par()$mar
  
  par(pty='s', mar=c(1,1,1,1))
  
  pal <- c(1,2)
  pal.alpha <- setNames(rgb(data.frame(t(col2rgb(pal))/255), alpha=.6),
                      names(pal))
  
  tcsdat <- suppressWarnings(tcs(dat))
    
  par(mar=c(1,1,1,1))
    
  test<-scatterplot3d(x=0, y=0, z=0, box=TRUE, pch=16, color=rgb(0,0,0,.5), 
                          xlim=c(-1.1,.3), ylim=c(-.25,.5), zlim=c(-.25,.65), 
                          axis=F, grid=F, angle=angle, scale.y=scale.y, 
                          cex.symbols=1, mar=c(2,1,1,1))
  
    #################################
    # INSERT POINT INFORMATION HERE #
    #################################
    
    test$points3d(tcsdat[1:(dim(dat)[1]/2), c('x', 'y', 'z')], 
                  pch=19, col=pal.alpha[1])
    test$points3d(tcsdat[(dim(dat)[1]/2+1):(dim(dat)[1]), c('x', 'y', 'z')], 
                  pch=19, col=pal.alpha[2])

    #########################
    # END POINT INFORMATION #
    #########################
    
    # Vertex coordinates
    uv<-test$xyz.convert(0,0,.75)
    blue<-test$xyz.convert((-1*sqrt(1.5)),(-1/(2*sqrt(2))),-.25)
    green<-test$xyz.convert(0,(1/sqrt(2)),-.25)
    red<-test$xyz.convert((.5*sqrt(1.5)),(-1/(2*sqrt(2))),-.25)
    no.uv<-test$xyz.convert(0,0,-.25)
    
    # Add text to vertices (can also be changed to points if preferred)
    #text(red$x,red$y,"l")
    #text(green$x,green$y,"m")
    #text(blue$x,blue$y,"s")
    #text(uv$x,uv$y,"u")
    
    points(uv$x, uv$y, pch=20, col="#984EA3", cex=1.5)
    points(red$x, red$y, pch=20, col="#E41A1C", cex=1.5)
    points(blue$x, blue$y, pch=20, col="#377EB8", cex=1.5)
    points(green$x, green$y, pch=20, col="#4DAF4A", cex=1.5)
    
    # Draw tetrahedron
    segcol <- grey(0.5)
    
    segments(uv$x,uv$y,red$x,red$y, col=segcol)
    segments(uv$x,uv$y,green$x,green$y, col=segcol)
    segments(uv$x,uv$y,blue$x,blue$y, col=segcol)
    segments(red$x,red$y,green$x,green$y, col=segcol)
    segments(green$x,green$y,blue$x,blue$y, col=segcol)
    segments(blue$x,blue$y,red$x,red$y, col=segcol)

    # Draw coordinates
    arrow1 <- 
      test$xyz.convert(
        matrix(c(
          0.4, 0.5, 0.4,
          0.5, 0.5, 0.4
        ), ncol=3, byrow=TRUE
        )
      )
    
    arrow2 <- 
      test$xyz.convert(
        matrix(c(
          0.4, 0.5, 0.4,
          0.4, 0.6, 0.4
        ), ncol=3, byrow=TRUE
        )
      )
    
    arrow3 <- 
      test$xyz.convert(
        matrix(c(
          0.4, 0.5, 0.4,
          0.4, 0.5, 0.5
        ), ncol=3, byrow=TRUE
        )
      )
    
    arrows(arrow1$x[1], arrow1$y[1], arrow1$x[2], arrow1$y[2], lwd=2, length=0.03)
    text(arrow1$x[2]*1.07,arrow1$y[2],"x")
    arrows(arrow2$x[1], arrow2$y[1], arrow2$x[2], arrow2$y[2], lwd=2, length=0.03)
    text(arrow2$x[2]*1.07,arrow2$y[2]*1.07,"y")
    arrows(arrow3$x[1], arrow3$y[1], arrow3$x[2], arrow3$y[2], lwd=2, length=0.03)
    text(arrow3$x[2],arrow3$y[2]*1.07,"z")
    
    cat('')
    #par(pty=oldpar.pty, mar=oldpar.mar)
}

# dichtcp(simdich())

