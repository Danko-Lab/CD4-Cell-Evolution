#
# Simple version of the transcript HMM.
#
library(rqhmm)
library(bigWig)
library(bitops)

make.data <- function(chrom, bwSet.plus, bwSet.minus, step = 50) {
  # for each strand:
  # 1. collect reads by step
  # 2. create combined position vector (position is used if it has reads in
  #    any of the tracks in the set)
  # 3. reverse is needed
  collect.data.set <- function(bwSet, is.minus) {
    values = lapply(bwSet, function(bw)
      abs(chromStepSum.bigWig(bw, chrom, step, defaultValue = 0)))
    # position mask
    position.mask = rep(FALSE, length(values[[1]])) ## Should all be the same length (same chromosome).
    for (reads in values)
      position.mask = position.mask | reads > 0
    positions = which(position.mask)
    n = length(positions)

    if (n == 0)
      return(NULL)

    if (!is.minus) {
      distances = positions[2:n] - positions[1:(n-1)]
      final.positions = positions[1:(n-1)]

      # make table
      tbl = rbind(final.positions, distances)
      missing = rep(0, length(final.positions))

      for (reads in values) {
        tbl = rbind(tbl, reads[final.positions])
        mrow = (reads[final.positions] == 0)*1
        missing = rbind(missing, mrow)
      }
      return(list(data = tbl, missing = missing))
    } else {
      distances = positions[2:n] - positions[1:(n-1)]

      final.positions = positions[2:n]

      # make table
      tbl = rbind(rev(final.positions), rev(distances))
      missing = rep(0, length(final.positions))
      
      for (reads in values) {
        tbl = rbind(tbl, rev(reads[final.positions]))
        mrow = (reads[final.positions] == 0)*1
        missing = rbind(missing, rev(mrow))
      }
      return(list(data = tbl, missing = missing))
    }
  }

  result = vector(mode="list", length=2)
  missing = vector(mode="list", length=2)
  res.plus = collect.data.set(bwSet.plus, is.minus=FALSE)
  result[[1]] = res.plus$data
  missing[[1]] = res.plus$missing

  res.minus = collect.data.set(bwSet.minus, is.minus=TRUE)
  result[[2]] = res.minus$data
  missing[[2]] = res.minus$missing

  return(list(data = result, missing = missing))
}



# HMM
makeHmm <- function() {

# Notes on Andre's code:
#
# data.shape: Sets the structure of the emission and covariate data passed to the HMM.  list(emission_dimenstions, covariates)
#   Emission_dimensions: Emissions of the model at each state.
#     * A vector of the dimenstions for each emission distribution.  
#     * >1 currently not supported in any meaningful way (distributions must be implemented).
#     * 2 sets emissions as a 2d distribution (e.g. a 2d guassian).
# 
#   Covariates: Additional data which can be applied to either emissions or transitions.  Each covariate slot can be assigned to multiple functions.
#      How these are used is up to the invdividual implementation of the transition/ emission functions.
#
# valid.transitions: Sets the transition functions availiable to the HMM.
#   NxN matrix (where N is the number of rows or columns in the matrix).
#   Transitions are represented by the position in the matrix.  Including a 0 prevents a transition during that pair of states.   
#
#   Valid transitions  must be numbered from 1 to N (where N is the number of valid positions in that row).  The number specifies the order in which 
#   transition.functions and values are specified.  Also the order of the transition vector in the set.transition.params.qhmm function.
# 
#	 More from Andre about these numbers: (Gchat on 10-6-2013).
#		"they are important in two ways
#		one is when setting parameters
#		depending on the transition function, these may be interpreted differently
#		in the case of the discrete transition function, they will just match the order in which you specify the transition probabilities
#		the second way in which these numbers are important, are when grouping states with respect to the transition functions
#		then, the transitions that have the same number are treated as being equivalent within the transition group"
#
# transition.fucntions
# emission.functions: For current support, see the bottom of: rqhmm/src/rqhmm.cpp.
#
# emission.groups: 
#  c(<slot>, <state1>, ..., <stateN>)
#   or
#  list(c(<slot1>, ..., <slotM>), c(<state1>, ..., <stateN>))

	n_states <- 32
	n_species <- 3

	#########################################################
	## Set up transition matrix.

	vtrans <- matrix(0, nrow=n_states, ncol=n_states)

	## Transitions from background to other states.
	vtrans[1, ] <- c(1, 3:8,2,rep(0,6), 10:11,9,0,0, 13:14,12,0,0, 16:17,15,0,0, 18, 19, 20)
	#                       ^--TTT      ^--2x species transcribed                ^--species specific gains.

	## Transitions from first set of shadow states (to TTT) (one species transcribed).
	vtrans[2, ] <- c(0, 1,0,0, 2,3,0, 4, rep(0,24))
	vtrans[3, ] <- c(0, 0,1,0, 2,0,3, 4, rep(0,24))
	vtrans[4, ] <- c(0, 0,0,1, 0,2,3, 4, rep(0,24))

	## Transitions from second set of shadow states (to TTT) (two species transcribed).
	vtrans[5, ] <- c(rep(0,4), 1,0,0, 2, rep(0,24))
	vtrans[6, ] <- c(rep(0,4), 0,1,0, 2, rep(0,24))
	vtrans[7, ] <- c(rep(0,4), 0,0,1, 2, rep(0,24))

	## Transitions from TTT back to normal.
	vtrans[8, ] <- c(2, rep(0,6), 1, 3:8, rep(0,18)) ## TTT
	#                ^--BBB          ^--back transition shadow states

	## Transitions from first set of back transition shadow states (two species transcribed).
	vtrans[9, ] <- c(4, rep(0,6), 0, 1,0,0, 2,3,0, rep(0,18))
	vtrans[10,] <- c(4, rep(0,6), 0, 0,1,0, 2,0,3, rep(0,18))
	vtrans[11,] <- c(4, rep(0,6), 0, 0,0,1, 0,2,3, rep(0,18))

	## Transitions from second set of back transition shadow states (one species transcribed).
	vtrans[12,] <- c(2, rep(0,10), 1,0,0, rep(0,18))
	vtrans[13,] <- c(2, rep(0,10), 0,1,0, rep(0,18))
	vtrans[14,] <- c(2, rep(0,10), 0,0,1, rep(0,18))

	## Transitions from shadow states (to BTT, TBT, and TTB)
	vtrans[15,] <- c(0, rep(0,13), 1,0, 2, rep(0,15))
	vtrans[16,] <- c(0, rep(0,13), 0,1, 2, rep(0,15))
	vtrans[17,] <- c(4, rep(0,13), 0,0, 1, 2,3, rep(0,13)) ## BTT
	vtrans[18,] <- c(2, rep(0,13), 0,0, 0, 1,0, rep(0,13))
	vtrans[19,] <- c(2, rep(0,13), 0,0, 0, 0,1, rep(0,13))

	vtrans[20,] <- c(0, rep(0,18), 1,0, 2, rep(0,10))
	vtrans[21,] <- c(0, rep(0,18), 0,1, 2, rep(0,10))
	vtrans[22,] <- c(4, rep(0,18), 0,0, 1, 2,3, rep(0,8)) ## TBT
	vtrans[23,] <- c(2, rep(0,18), 0,0, 0, 1,0, rep(0,8))
	vtrans[24,] <- c(2, rep(0,18), 0,0, 0, 0,1, rep(0,8))

	vtrans[25,] <- c(0, rep(0,23), 1,0, 2, rep(0,5))
	vtrans[26,] <- c(0, rep(0,23), 0,1, 2, rep(0,5))
	vtrans[27,] <- c(4, rep(0,23), 0,0, 1, 2,3, rep(0,3)) ## TTB
	vtrans[28,] <- c(2, rep(0,23), 0,0, 0, 1,0, rep(0,3))
	vtrans[29,] <- c(2, rep(0,23), 0,0, 0, 0,1, rep(0,3))

	## Transitions to species-specific gains.
	vtrans[30,] <- c(2, rep(0,28), 1, 0, 0) ## TBB
	vtrans[31,] <- c(2, rep(0,28), 0, 1, 0) ## BTB
	vtrans[32,] <- c(2, rep(0,28), 0, 0, 1) ## BBT


	#####################################################
	## Create an emissions group matrix.
	
	slot1_bg <- c(1,3,4,7,11,13,14,15,16,17,18,19,21,24,26,29,31,32)
	slot1_tr <- c(1:32)[!(c(1:32) %in% slot1_bg)]

	slot2_bg <- c(1,2,4,6,10,12,14,15,18,20,21,22,23,24,25,28,30,32)
	slot2_tr <- c(1:32)[!(c(1:32) %in% slot2_bg)]

	slot3_bg <- c(1,2,3,5, 9,12,13,16,19,20,23,25,26,27,28,29,30,31)
	slot3_tr <- c(1:32)[!(c(1:32) %in% slot3_bg)]

	emissions <- new.emission.groups(n_states, 1+n_species)
	emissions <- add.emission.groups(emissions, states= c(1:n_states), slots= rep(1, n_states))
	emissions <- add.emission.groups(emissions, states= c(slot1_bg,slot2_bg,slot3_bg),
								slots= c(rep(2, NROW(slot1_bg)), rep(3, NROW(slot2_bg)), rep(4, NROW(slot3_bg))))
	emissions <- add.emission.groups(emissions, states= c(slot1_tr,slot2_tr,slot3_tr),
								slots= c(rep(2, NROW(slot1_tr)), rep(3, NROW(slot2_tr)), rep(4, NROW(slot3_tr))))
	emissions.functions <- c("geometric", rep("neg_binomial", n_species))

	#####################################################
	## Create the HMM.
	
	hmm = new.qhmm(data.shape= list(c(1, rep(1, n_species)), NULL), support.missing = TRUE,
	  valid.transitions= vtrans,
	  transition.functions= rep("discrete", n_states),
	  emission.functions= as.list(rep(list(emissions.functions), n_states)), emission.groups= emissions)

	#####################################################
	## Set up model parameters.

	g <- 1e-3 ## \gamma --> Parameter for a complete gain or loss.

	b <- 1e-2 ## \beta  --> Parameter for a transition from B->T.
	a <- 1e-2 ## \alpha --> Parameter for a change in start or end site.
	d <- 1e-2 ## \delta --> Parameter specifying the 'length' of 'extension' in start site.

	e <- 1e-2 ## \eta --> Transition from T->B.
	ap<- a    ## \alpha^prime --> Parameter for a change in end site.
	dp<- d    ## \delta --> Parameter specifying the 'length' of 'extension' in end site.


	## Phylogenetic tree.
	H <- 6/(6+6+25)  ## branch length for human/ (h+c+m)
	C <- 6/(6+6+25)  ## branch length for chimp/ (h+c+m)
	M <- 25/(6+6+25) ## branch length for Rmac./ (h+c+m)

	## Set the initial transition parameters...
	set.transition.params.qhmm(hmm, 1, c(1-b, 
			b*(1-a-g-a*g), a*H*b/2, a*C*b/2, a*M*b/2, a*M*b/2, a*C*b/2, a*H*b/2, #TTT
			g*H*b/2, g*a*(M+C)*b/4, g*a*(C+M)*b/4,                               #BTT
			g*C*b/2, g*a*(M+C)*b/4, g*a*(C+M)*b/4,                               #TBT
			g*M*b/2, g*a*(H+C)*b/4, g*a*(C+H)*b/4,                               #TTB
			g*H*b/2, g*C*b/2, g*M*b/2                             ),  fixed=rep(TRUE, 20)) #TBB,BTB,BBT
			
	set.transition.params.qhmm(hmm, 8, c( 1-e, e*(1-ap),
			e*ap*M/2, e*ap*C/2, e*ap*H/2, e*ap*H/2, e*ap*C/2, e*ap*M/2 ),  fixed=rep(TRUE, 8))
			
	## Intermediate states with T in 1 species -> TTT|(T in 2 species).
	  set.transition.params.qhmm(hmm, c(2:4), c(1-d, d*a/2, d*a/2, d*(1-a)), fixed=rep(TRUE, 4))
	  
#	## Intermediate states with T in 1 species -> TTT|(other knockout state).
	  set.transition.params.qhmm(hmm, c(5:7,15:16,20:21,25:26), c(1-d, d), fixed=rep(TRUE, 2))

	## Intermediate states with T in 2 species -> BBB|(T in 1 species).
	  set.transition.params.qhmm(hmm, c(9:11), c(1-dp, dp*ap/2, dp*ap/2, dp*(1-ap)), fixed=rep(TRUE, 4))

	## Intermediate states with T in 1 species -> BBB.
	  set.transition.params.qhmm(hmm, c(12:14,18:19,23:24,28:29), c(1-dp, dp), fixed=rep(TRUE, 2))

	## For species-specific losses :OR: gains on HC branch.
	  set.transition.params.qhmm(hmm, c(17,22,27), c(1-e, e*ap/2, e*ap/2, e*(1-ap)), fixed=rep(TRUE, 4))

	## For species-specific gains on terminal branches.
	  set.transition.params.qhmm(hmm, c(30:32), c(1-e,e), fixed=rep(TRUE, 2))

	#################################
	## Set initial probabilities.
	set.initial.probs.qhmm(hmm, c(1, rep(0, n_states - 1))) # start with background

	#################################
	## Set emissions params.
	set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
    set.emission.params.qhmm(hmm, c(2:n_states), 1/20, slot = 1) ## All other states have at least one sequence transcribed.
	
	back_emissions <- c(5, 1/5) ## NOTE: Setting a third parameter equal to the 'mean' fixes the mean.
	tran_emissions <- c(10, 3/10)

   set.emission.params.qhmm(hmm, slot1_bg, back_emissions, slot = 2)
   set.emission.params.qhmm(hmm, slot1_tr, tran_emissions, slot = 2)
   set.emission.params.qhmm(hmm, slot2_bg, back_emissions, slot = 3)
   set.emission.params.qhmm(hmm, slot2_tr, tran_emissions, slot = 3)
   set.emission.params.qhmm(hmm, slot3_bg, back_emissions, slot = 4)
   set.emission.params.qhmm(hmm, slot3_tr, tran_emissions, slot = 4)

   return(hmm)
}

#########
#set.transition.option.qhmm ## Set the tree...

decode_hmm <- function(chrom, hmm, full.data, train.data, missing.lst, step = 50, min.len = NA) {

  # common info
  n.tracks = dim(train.data[[1]])[1] - 1
  next.idx = 1

  print.pattern.stats <- function(stats) {
    cat.pattern <- function(pat) {
      for (i in 1:n.tracks) {
        if (bitAnd(pat, 2^(i - 1)) > 0)
          cat("#")
        else
          cat(" ")
      }
    }

    m = 2^n.tracks
    for (k in 1:m) {
      cat.pattern(k)
      cat(" :", stats[k], "\n")
    }
  }
  
  decode.path <- function(path, positions, strand) {
    tcount.stats = rep(0, n.tracks)
    pattern.stats = rep(0, 2^n.tracks)
    
    blocks = path.blocks.qhmm(path, 2:(hmm$n.states))

    cat("# blocks:", dim(blocks)[2], "\n")

    coords = rbind(positions[blocks[1,]], positions[blocks[2,]])
    if (strand == '-')
      coords = rbind(positions[blocks[2,]], positions[blocks[1,]])

    len = coords[2,] - coords[1,] + 1
    if (!is.na(min.len)) {
      coords = coords[, len > min.len]
      blocks = blocks[, len > min.len]
    }

#    coords = t(coords)*step
#    coords[,1] = coords[,1] - 1 # convert to zero-based, right-open
    N = dim(blocks)[2]

    cat("# blocks:", N, "\n")
    
    #
    # now need to produce the sub-tracks by examining the state path
    # within the block
    starts = NULL
    ends = NULL
    names = NULL
    for (i in 1:N) {
      idx.start = blocks[1, i]
      idx.end = blocks[2, i]

      subpath = path[idx.start:idx.end]

      tcount = 0
      pattern = 0
      for (j in 1) { #:n.tracks) {

        edge.offset = length(subpath) #1
        start.ij = positions[idx.start]
        end.ij = positions[idx.start + edge.offset - 1]
        if (strand == '-') {
          start.ij = positions[idx.start + edge.offset - 1]
          end.ij = positions[idx.start]
        }

        start.ij.coord = start.ij * step - 1
		end.ij.coord = (end.ij + 1) * step - 1
        if (is.na(min.len) || (end.ij.coord - start.ij.coord) > min.len) {
          starts = c(starts, start.ij.coord)
          ends = c(ends, end.ij.coord)
          names = c(names, paste("T:", chrom, ":", next.idx, ":", j, sep=''))
          tcount = tcount + 1
          pattern = pattern #+ bmask
        }
      }
      tcount.stats[tcount] = tcount.stats[tcount] + 1
      pattern.stats[pattern] = pattern.stats[pattern] + 1
      if (tcount > 0)
        next.idx <<- next.idx + 1
    }

    print(tcount.stats)
    print.pattern.stats(pattern.stats)
    
    return(data.frame(chrom = chrom, start = starts, end = ends, name = names, score = 0, strand = strand))

    # return number of groups
  }

  path.plus = viterbi.qhmm(hmm, train.data[[1]], missing = missing.lst[[1]])
  path.minus = viterbi.qhmm(hmm, train.data[[2]], missing = missing.lst[[2]])
  
  bed.plus = decode.path(path.plus, full.data[[1]][1,], "+")
  bed.minus = decode.path(path.minus, full.data[[2]][1,], "-")

  print(next.idx - 1)
  
  return(rbind(bed.plus, bed.minus))
}

# svn switch --relocate $SVN_ROOT/phast/trunk svn+ssh://cgd24@compgen.bscb.cornell.edu/usr/local/svnroot/phast/trunk .
 # Data: H1, C3, R2.
 Hp <- load.bigWig("/usr/projects/GROseq/CD4/Alignments/CD4-U_plus.bw")
 Hm <- load.bigWig("/usr/projects/GROseq/CD4/Alignments/CD4-U_minus.bw")
 Cp <- load.bigWig("/usr/projects/GROseq/NHP/Alignments_3rdPrep/C4-U.bed.gz_plus.hg19.bw")
 Cm <- load.bigWig("/usr/projects/GROseq/NHP/Alignments_3rdPrep/C4-U.bed.gz_minus.hg19.bw")
 Mp <- load.bigWig("/usr/projects/GROseq/NHP/Alignments_3rdPrep/M3-U.bed.gz_plus.hg19.bw")
 Mm <- load.bigWig("/usr/projects/GROseq/NHP/Alignments_3rdPrep/M3-U.bed.gz_minus.hg19.bw")

 all.data = make.data("chr22", list(Hp, Cp, Mp), list(Hm, Cm, Mm))

 ## SANITY CHECK
 #head(t(all.data$data[1][[1]])) ## SANITY CHECK
 
 dataset.train = lapply(all.data$data, function(data) { data[2:5,] })
 hmm <- makeHmm()
 em.qhmm(hmm, dataset.train, missing.lst = all.data$missing, n_threads = 4)
 collect.params.qhmm(hmm)
 
 bed <- decode_hmm("chr22", hmm, all.data$data, dataset.train, all.data$missing, step = 50, min.len = NA)
 write.table(bed, "tmp.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
 
 path_plus = viterbi.qhmm(hmm, all.data$data[[1]][c(2:5),], missing = all.data$missing[[1]])
 path_minus = viterbi.qhmm(hmm, all.data$data[[2]][c(2:5),], missing = all.data$missing[[2]])

 ## SANITY CHECK
 head(cbind( (t(all.data$data[1][[1]])), path_plus ), 100) ## SANITY CHECK
 
# posterior decoding
fw = forward.qhmm(hmm, rolls)
bk = backward.qhmm(hmm, rolls)

logPx = attributes(fw)$loglik

posterior = exp(fw + bk - logPx)

