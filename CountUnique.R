##
##  CountUnique.R  --  File contains a few conveninece functions that help counting, and 
##                     identifying unique transcripts.
##
##  4/25/10 -- Should be all set!  They agree with one another, and pass all of the DUH tests.
##

require()

## ... 
isIn <- function(f, p) {
  
}


## Counts the number of unique transcripts in a mixture, using the GROseq packages' AssociateWithInterval.
countUnique <- function(UP) {
 UPFirstUnique <- UP
 NUM <- 1
 while( NROW(UPFirstUnique) > 1 ) {
   ## Each time, this removes index=1, and any overlaps.
   UPFirstUnique <- UPFirstUnique[is.na(associateWithInterval(UPFirstUnique[1,], UPFirstUnique)),]
   NUM<- NUM+1  ## Increment each time it removes.
 }
 print(NUM)
}

## Returns a vector, the same size as G, which represents whether each G is the first of a unique genomic region.
is.firstUnique <- function(G) {
  fu <- rep(FALSE, NROW(G))
  fu[1] <- TRUE
  for(i in 2:NROW(G)) {
    fu[i] <- (sum(!is.na(associateWithInterval(f=G[i,], p=G[c(1:(i-1)),]))) == 0)
  }
  return(fu)
}


