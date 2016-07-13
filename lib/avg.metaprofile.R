## Alternative to Andre's bigWig metaProfile function.
## only change is that op= "avg" insteand of "sum".
##
## Provides more interpretable Y-axis values for some data types (e.g., phylo-P).

#
# Convinience function
#

# Profile object structure:
# . name
# . x0 (in units of number of steps)
# . step (size in bp)
# . top, middle, bottom vectors
#

avg.metaprofile.bigWig <- function(bed, bw.plus, bw.minus = NULL, step = 1, name = "Signal", matrix.op = NULL, profile.op = subsampled.quantiles.metaprofile, ...) {
  #
  # 1. collect data
  N = dim(bed)[2]
  mat = NULL
  if (N >= 6) {
    stopifnot(bw.minus != NULL)
    mat = bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed, step, op = "avg", abs.value = TRUE, as.matrix = TRUE, follow.strand = TRUE)
  } else {
    if (!is.null(bw.minus))
      warning("bw.minus != NULL but BED contains no strand information")
    mat = bed.step.bpQuery.bigWig(bw.plus, bed, step, op = "avg", as.matrix = TRUE)
  }
  
  #
  # 2. apply matrix transformation
  if (!is.null(matrix.op))
    mat = matrix.op(mat, ...)
  
  #
  # 3. apply profile operation
  result = profile.op(mat, ...)
  
  #
  # 4. create result
  X0 = 0 # can't really tell what X0 was from the input arguments
  
  res = c(list(name = name, X0 = X0), result)
  
  class(res) <- "metaprofile"
  return(res)
}

