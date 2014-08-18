source("polyA.vs.SS5.qhmm.R", chdir=TRUE)
hg19 = twobit.load("/gbdb/hg19/hg19.2bit")
rheMac3 = twobit.load("/gbdb/rheMac3/rheMac3.2bit")
panTro4 = twobit.load("/gbdb/panTro4/panTro4.2bit")
scan.length = 1500 # 1.5 kbp

collect.sequence <- function(twobit, bed, seq.length = 1500) {
  N = dim(bed)[1]
  result = vector(mode="list", length = N)

  is.minus = bed[,6] == '-'
  start = bed[,2]
  end = bed[,2] + seq.length
  start[is.minus] = bed[is.minus, 3] - seq.length
  end[is.minus] = bed[is.minus, 3]

  bed.clipped = data.frame(bed[,1], start, end, bed[,4:6])
  
  foreach.bed(bed.clipped, function(i, chrom, start, end, strand) {
    seq = twobit.sequence(twobit, chrom, start, end)
      
      if (strand == '-')
        seq = twobit.reverse.complement(seq)

      result[[i]] <<- twobit.sequence.to.integer(seq)
  })

  return(result)
}

path.probs.summary <- function(hmm) {
  # get re-normalized outgoing probabilities
  path.probs = get.transition.params.qhmm(hmm, 1)[2:6]
  path.probs = path.probs / sum(path.probs)

  # combine and label paths
  probs = c(
    path.probs[1],
    path.probs[2] + path.probs[3],
    path.probs[4] + path.probs[5])

  names(probs) <- c("None", "SS5 first", "pA first")
  return(probs)
}


#
# predict stability based on posterior probability of being in
# the pA first paths
#

unstable.score <- function(hmm, data, n.threads = 1) {
  sapply(data, function(seq) {
    post = posterior.from.state.qhmm(hmm, 1, seq, n_threads = n.threads)

    # outgoing transitions: B, End, SS5.A1, SS5, pA, pA.A2
    out.probs = rowSums(post)[2:6] # drop B
    out.probs = out.probs / sum(out.probs) # renormalize

    # unstable prob := prob of finding a pA first
    sum(out.probs[4:5])
  })
}

