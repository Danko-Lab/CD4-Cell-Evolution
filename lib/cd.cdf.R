## Simple R script to draw a CDF using base graphics.
## Point is to make CDFs smaller than those created by 
## R's internal function.  Can import into Illustrator/ Canvas.
## 

cd.cdf <- function(x, npoints=200, add=FALSE, ...) {
  x.axis <- unique(sort(x))
  y.axis <- sapply(x.axis, function(i) {sum(i>=x, na.rm=TRUE)/sum(!is.na(x))})

  if((NROW(x.axis) > (2*npoints)) & !is.na(npoints)) {
    indx <- c(1, which((2:(NROW(x.axis)-1) %% round(NROW(x.axis)/npoints)) == 0), NROW(x.axis))
    x.axis <- x.axis[indx]
    y.axis <- y.axis[indx]
  }
  
  if(add==TRUE) {
	points(y.axis~x.axis, ...)
  }
  else {
	plot(y.axis~x.axis, ylim=c(0,1), ...)
  }
}

#cd.cdf(rnorm(5000), col="dark green", lwd=3, type="l"); abline(h=0.5); abline(v=0)
#cd.cdf(rnorm(5000)+0.5, add=TRUE, col="dark red", lwd=3, type="l"); 
