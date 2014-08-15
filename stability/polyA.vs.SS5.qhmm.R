#
# reimplementation of the polyA vs SS5 HMM using the QHMM package
#
library(grocap)
library(rqhmm)

#
# B -> :: 
#  : End
#  : SS5 -> sink
#  : SS5 -> B -> pA -> sink
#  : pA -> sink
#  : pA -> B -> SS5 -> sink
# sink(B) -> End
#
# refactored:
# B -> ::
#  : End
#  A1: SS5 -> sink
#  A2: pA -> sink
#  : SS5 -> B(2) -> (A2)
#  : pA -> B(3) -> (A1)
# sink(B) -> End
#
make.ss5.pA.qhmm <- function(back, ss5, polyA) {
  L1 = dim(ss5)[2]
  L2 = dim(polyA)[2]

  N = 5 + 2*L1 + 2*L2 # B, End, sink(B), B(2), B(3) + SS5, SS5(A1), pA, pA(A2)

  # state numbers
  B = 1
  End = 2
  SS5.A1 = 3
  SS5 = 3 + L1
  B2 = SS5 + L1
  pA = B2 + 1
  B3 = pA + L2
  pA.A2 = B3 + 1
  sink = pA.A2 + L2

  # make transition map
  tm = matrix(data=0, nrow=N, ncol=N)
  tm[B, c(B, End, SS5.A1, SS5, pA, pA.A2)] = 1:6
  tm[End, End] = 1 # just make the HMM code happy

  for (i in 1:(L1-1))
    tm[SS5.A1 + i - 1, SS5.A1 + i] = 1
  tm[SS5.A1 + L1 - 1, sink] = 1

  for (i in 1:(L1-1))
    tm[SS5 + i - 1, SS5 + i] = 1
  tm[SS5 + L1 - 1, B2] = 1

  tm[B2, c(B2, pA.A2)] = 1:2

  for (i in 1:(L2-1))
    tm[pA + i - 1, pA + i] = 1
  tm[pA + L2 - 1, B3] = 1

  tm[B3, c(B3, SS5.A1)] = 1:2

  for (i in 1:(L2-1))
    tm[pA.A2 + i - 1, pA.A2 + i] = 1
  tm[pA.A2 + L2 - 1, sink] = 1

  tm[sink, c(sink, End)] = 1:2

  #
  hmm <- new.qhmm(list(1, NULL),
                  tm, rep("discrete", N),
                  as.list(rep("discrete", N)))

  #return(hmm)
  
  # set transition parameters
  # NOTE: only need to set those that are different from 1
  #       but I want to make parameters fixed, so ...
  set.transition.params.qhmm(hmm, B, c(0.997, rep(0.003 / 5, 5))) # not fixed
  set.transition.params.qhmm(hmm, B2, c(0.997, 0.003), fixed = c(T, T))
  set.transition.params.qhmm(hmm, B3, c(0.997, 0.003), fixed = c(T, T))
  set.transition.params.qhmm(hmm, sink, c(0.997, 0.003), fixed = c(T, T))

  other.states = setdiff(1:N, c(B, B2, B3, sink))
  set.transition.params.qhmm(hmm, other.states, 1, fixed = T)

  # set emission parameters
  # background like states
  set.emission.params.qhmm(hmm, c(B, B2, B3, sink), c(0, back), fixed = rep(T, 5))
  # end state
  set.emission.params.qhmm(hmm, End, c(1, rep(0, 4)), fixed = rep(T, 5))
  # SS5 1 motif
  for (i in 1:L1)
    set.emission.params.qhmm(hmm, SS5.A1 + i - 1, c(0, ss5[,i]), fixed = rep(T, 5))
  # SS5 2 motif
  for (i in 1:L1)
    set.emission.params.qhmm(hmm, SS5 + i - 1, c(0, ss5[,i]), fixed = rep(T, 5))
  # pA
  for (i in 1:L2)
    set.emission.params.qhmm(hmm, pA + i - 1, c(0, polyA[,i]), fixed = rep(T, 5))
  # pA.A2
  for (i in 1:L2)
    set.emission.params.qhmm(hmm, pA.A2 + i - 1, c(0, polyA[,i]), fixed = rep(T, 5))

  # set initial state probs
  set.initial.probs.qhmm(hmm, c(1, rep(0, N - 1)))
  
  return(hmm)
}

#
# PWMs
#

load("../../splicing/splicing.pwms.Rdata")

# poly A
polyA.patterns = list(
  list(58.2, c(1, 1, 4, 1, 1, 1)),
  list(14.9, c(1, 4, 4, 1, 1, 1)),
  list( 2.7, c(1, 3, 4, 1, 1, 1)),
  list( 3.2, c(4, 1, 4, 1, 1, 1)),
  list( 1.3, c(2, 1, 4, 1, 1, 1)),
  list( 1.3, c(3, 1, 4, 1, 1, 1)),
  list( 1.7, c(1, 1, 4, 1, 4, 1)),
  list( 1.2, c(1, 1, 4, 1, 2, 1)),
  list( 0.7, c(1, 1, 4, 1, 3, 1)),
  list( 0.8, c(1, 1, 1, 1, 1, 3)),
  list( 0.6, c(1, 2, 4, 1, 1, 1)))

polyA.w = sapply(polyA.patterns, function(lst) lst[[1]])

polyA.col <- function(col) {
  bases = sapply(polyA.patterns, function(lst) lst[[2]][col])
  
  sapply(1:4, function(b) sum(polyA.w[bases == b]) / sum(polyA.w))
}

polyA.pwm = sapply(1:6, polyA.col)
ss5.pwm = res.genes$ss5

#
# data
make.ss5.pA.hmm.data <- function(seqlst) {
  lapply(seqlst, function(seq) c(seq + 1, 1))
}


res = make.ss5.pA.qhmm(rep(1/4, 4), ss5.pwm, polyA.pwm)
