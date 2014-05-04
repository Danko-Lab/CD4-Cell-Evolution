#
# Normalizes two samples by creating an empirical distribution
# and subsampling.
#
# Returns indices for each sample.
#

norm.subsample <- function(a, b, nBins=100, nsamp=10000, boot.replace=TRUE, plot.cdf=FALSE) {

 ## Choose one of the datasets, DS.
 ## Descritise DS.  Get normalized numbers in each bin.
 rs <- range(a, b)
 binInc  <- (rs[2]-rs[1])/nBins
 binMids <- seq(rs[1]-binInc,rs[2]+binInc,binInc)

 ## Get a combined distribution for a and b.
 binFreq <- sapply(binMids, function(x) { (sum(a<(x+binInc/2) & a>(x-binInc/2))+sum(b<(x+binInc/2) & b>(x-binInc/2))) /	(NROW(a)+NROW(b)) })
 
 ## Compute probabilities for each element in a and b.
 a.sweight <- rep(0,NROW(a))
 b.sweight <- rep(0,NROW(b))
 for(i in c(1:NROW(binMids))) {
   bM <- binMids[i]
   incl <- a<(bM+binInc/2) & a>(bM-binInc/2)
   a.sweight[incl] <- binFreq[i]/sum(incl) ## Distribute frequency of that bin among ## in the dataset...
 
   incl <- b<(bM+binInc/2) & b>(bM-binInc/2)
   b.sweight[incl] <- binFreq[i]/sum(incl) ## Distribute frequency of that bin among ## in the dataset...
 }

 ## Subsample with a probability proportional to the combination of a and b.
 a.sIndx <- sample(c(1:NROW(a)), nsamp, prob= a.sweight, replace=boot.replace)
 b.sIndx <- sample(c(1:NROW(b)), nsamp, prob= b.sweight, replace=boot.replace)

 ## Sanity check that subsamplin worked as advertised...
 if(plot.cdf) {
   par(mfrow=c(1,1))
   plot(ecdf(a), xlim=rs, col="pink", main="CDF before (light) and after (dark) subsampling.")
   plot(ecdf(b), col="light blue", add=TRUE)
   plot(ecdf(a[a.sIndx]), col="red", add=TRUE)
   plot(ecdf(b[b.sIndx]), col="blue", add=TRUE)
 }

 list(s1= a.sIndx, s2= b.sIndx)

}

## Generalizes norm.subsample for an arbitrary number of classes.
## l is a list of vectors.  Each represents a distinct class.
norm.subsample.n <- function(l, nBins=100, nsamp=1000, boot.replace=TRUE, plot.cdf=FALSE) {

 ## Descritise DS.
 rs <- range(unlist(l))
 binInc  <- (rs[2]-rs[1])/nBins
 binMids <- seq(rs[1]-binInc,rs[2]+binInc,binInc)

 ## Get a combined distribution over all of the classes.
 ne <- NROW(unlist(l))
 binFreq <- sapply(binMids, function(x) { (sum( sapply(l, function(y) {sum(y<(x+binInc/2) & y>(x-binInc/2))}) )/ ne)  })
 
 ## Compute the probabilities for each element in each list.
 sweight <- l
 for(i in 1:NROW(binMids)) {
   bM <- binMids[i]
   for(j in 1:NROW(l)) {
     incl <- l[[j]]<(bM+binInc/2) & l[[j]]>(bM-binInc/2)
     sweight[[j]][incl] <- binFreq[i]/sum(incl)
   }
 }

 ## Subsample with a probability proposrtional to the combination of all elements of l.
 sIndx <- l
 for(j in 1:NROW(l)) {
  sIndx[[j]] <- sample(c(1:NROW(l[[j]])), nsamp, prob=sweight[[l]], replace=boot.replace)
 }

 ## Sanity check that subsampling worked ...
 if(plot.cdf) {
   plot(ecdf(l[[1]]), xlim=rs, main="CDF before and after subsampling", col="gray")
   plot(ecdf(l[[1]][sIndx[1]]), add=TRUE, col="black")
   for(j in 2:NROW(l)) {
     plot(ecdf(l[[j]]), col="gray", add=TRUE)
     plot(ecdf(l[[j]][sIndx[j]]), add=TRUE, col="black")
   }
 }

 return(sIndx)
}

