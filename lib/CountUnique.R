##
##  CountUnique.R  --  File contains a few conveninece functions that help counting, and 
##                     identifying unique transcripts.
##
##  4/25/10 -- Should be all set!  They agree with one another, and pass all of the DUH tests.
##

## ... 
isIn <- function(f, p) {
  stopifnot(NROW(f) == 1)
  chrom<- f[1,1]
  chromStart <- f[1,2]
  chromEnd   <- f[1,3]
  strand<- f[1,6]
  indx <- as.character(p[,1]) == chrom & p[,6] == strand
  sum(p[indx,3] > chromStart & p[indx,2] < chromEnd)
}


## Counts the number of unique transcripts in a mixture, using the GROseq packages' AssociateWithInterval.
countUnique <- function(UP) {
  return(sum(is.firstUniue(UP)))
}

## Returns a vector, the same size as G, which represents whether each G is the first of a unique genomic region.
is.firstUnique <- function(G) {
  fu <- rep(FALSE, NROW(G))
  fu[1] <- TRUE
  for(i in 2:NROW(G)) {
    fu[i] <- isIn(G[i,], G[c(1:(i-1)),]) == 0
  }
  return(fu)
}


