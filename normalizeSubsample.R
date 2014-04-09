#
# Normalizes two samples by creating an empirical distribution
# and subsampling.
#
# Returns indices for each sample.
#

test <- function() {
 a <- log(body$mean_Human)
 b <- log(pause$mean_Human)

 par(mfrow=c(3,1))
 hist(a, 50)
 hist(b, 50)
 plot(binMids, binFreq, type="l")

}

norm.subsample <- function(a, b, nBins=100, nsamp=10000, boot.replace=TRUE, plot.cdf=FALSE) {

 ## Choose one of the datasets, DS.
 ## Descritise DS.  Get normalized numbers in each bin.
 rs <- range(a, b)
 binInc  <- (rs[2]-rs[1])/nBins
 binMids <- seq(rs[1]-binInc,rs[2]+binInc,binInc)
 binFreq <- sapply(binMids, function(x) { (sum(a<(x+binInc/2) & a>(x-binInc/2))+sum(b<(x+binInc/2) & b>(x-binInc/2))) /	(NROW(a)+NROW(b)) })
 
## SANITY CHECK
# par(mfrow=c(3,1))
# hist(a, 50); points(binMids, binFreq*NROW(a), type="l", col="dark gray")
# hist(b, 50); points(binMids, binFreq*NROW(b), type="l", col="dark gray")
# plot(binMids, binFreq, type="l")

  a.sweight <- rep(0,NROW(a))
  b.sweight <- rep(0,NROW(b))
  for(i in c(1:NROW(binMids))) {
     bM <- binMids[i]
     incl <- a<(bM+binInc/2) & a>(bM-binInc/2)
     a.sweight[incl] <- binFreq[i]/sum(incl) ## Distribute frequency of that bin among ## in the dataset...
	 
     incl <- b<(bM+binInc/2) & b>(bM-binInc/2)
     b.sweight[incl] <- binFreq[i]/sum(incl) ## Distribute frequency of that bin among ## in the dataset...
  }
  a.sIndx <- sample(c(1:NROW(a)), nsamp, prob= a.sweight, replace=boot.replace)
  b.sIndx <- sample(c(1:NROW(b)), nsamp, prob= b.sweight, replace=boot.replace)

## SANTIY CHECK 2
# par(mfrow=c(3,1))
# hist(a[a.sIndx], 50); points(binMids, binFreq*nsamp, type="l", col="dark gray")
# hist(b[b.sIndx], 50); points(binMids, binFreq*nsamp, type="l", col="dark gray")
# plot(binMids, binFreq, type="l")

## Sanity check that subsamplin worked as advertised...
 if(plot.cdf) {
   par(mfrow=c(1,1))
   plot(ecdf(a), xlim=rs, col="pink", main="CDF before (light) and after (dark) subsampling.")
   plot(ecdf(b), col="light blue", add=TRUE)
   plot(ecdf(a[a.sIndx]), col="red", add=TRUE)
   plot(ecdf(b[b.sIndx]), col="blue", add=TRUE)
 }
}
