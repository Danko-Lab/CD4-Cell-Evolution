## 
## Centers dREG sites on the positions of highest read density.

center_dREG_site <- function(bed, bw_plus_file, bw_minus_file, readThresh=10) {
  require(bigWig)

#  readThresh <- 10 ## Mainly want to eliminate passing position bed[1,i] for regions without any read density.
  size=100
  overlap=25
  nwindows <- size/overlap

  ## Process command arguments.
#  bed_file <- tss#args[1]
#  bw_plus_file <- "../AllData/All_Merge/H-U_plus.bw"#args[3]
#  bw_minus_file <- "../AllData/All_Merge/H-U_minus.bw"#args[4]

  ## Load the bed and bigWigs. 
  bw_plus <- load.bigWig(bw_plus_file)
  bw_minus<- load.bigWig(bw_minus_file)

  ## Get the start position of 100bp window with the max counts inside each dREG site.
  chrom <- character()
  chromStart <- integer()
  chromEnd   <- integer()

  for(i in 1:NROW(bed)) {
    pad <- overlap- ((bed[i,3]-bed[i,2]) %% overlap) ## Padding, to make the window a multiple of the step size.
    bases  <- seq(bed[i,2], bed[i,3]+pad, overlap)

    countsp <- abs(step.bpQuery.bigWig(bw_plus, chrom=as.character(bed[i,1]), start= bed[i,2], end= bed[i,3]+pad, step=overlap))
    countsp_SUM <- sapply(1:(length(countsp)-nwindows), function(x) {sum(countsp[x:(x+nwindows)])})
    stopifnot(NROW(bases)-1-nwindows == NROW(countsp_SUM)) ## SANTIY CHECK

    countsm <- abs(step.bpQuery.bigWig(bw_minus, chrom=as.character(bed[i,1]), start= bed[i,2], end= bed[i,3]+pad, step=overlap))
    countsm_SUM <- sapply(1:(length(countsm)-nwindows), function(x) {sum(countsm[x:(x+nwindows)])})
    stopifnot(NROW(bases)-1-nwindows == NROW(countsm_SUM)) ## SANTIY CHECK

    if(max(countsp_SUM+countsm_SUM) > readThresh) {
      chromStart  <- c(chromStart, bases[1]+(which.max(countsp_SUM+countsm_SUM)-1-nwindows)*overlap)
      chromEnd    <- c(chromEnd, bases[1]+(which.max(countsp_SUM+countsm_SUM)+nwindows)*overlap)
      chrom       <- c(chrom, as.character(bed[i,1]))
      # paste(chrom, paste(chromStart, chromEnd, sep="-"), sep=":")
    } else { ## Take center of site.
      center <- bed[i,2]+round((bed[i,3]-bed[i,2])/2)
      chromStart <- c(chromStart, center-size)
      chromEnd   <- c(chromEnd,   center+size)
      chrom       <- c(chrom, as.character(bed[i,1]))
    }
  
  }

  data <- rbind(data.frame(chrom= chrom, start= chromStart, end= chromEnd, name= ".", score= "0", strand= "."))
  return(data)
}
