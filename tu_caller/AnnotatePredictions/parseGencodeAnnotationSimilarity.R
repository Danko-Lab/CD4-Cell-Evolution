args <- commandArgs(trailingOnly = TRUE)

## Command line arguments, passed with --args
matchfile <- args[1] #"TMP/GENCODE.overlap.stats"
transfile <- args[2] #"../bedCG/MakeFinalVersion/FinalTranscripts.bed"
annotfile <- args[3] #"GENCODE/gencode.v7.transcript.l_1_2_3.tsv"
outfileMapped <- args[4] #"TMP/GENCODE.Overlap.bed"
outfileUnmapped <- args[5] #"TMP/GENCODE.NoOverlap.bed"

## Read transcript overlap stats ...
gcSim <- read.table(matchfile, header=F)

## Read original data table with annotation data and transcript data.
GENCODE <- read.table(annotfile)
TRANSCRIPT <- read.table(transfile)

## Find the (semi)-unique ID of the GENCODE transcript with the best match to the HMM detected transcript.
#bestAnnotation <- unlist(lapply(unique(gcSim$V1), function(x) {
bestAnnotation <- NULL
for(x in unique(gcSim$V1)) {
        indx <- which(gcSim$V1 == x)
	ofIndx <- which.max(gcSim$V6[indx])
	bestAnnotation <- c(bestAnnotation, as.character(gcSim$V2[indx][ofIndx]))
	#return(gcSim$V2[indx][ofIndx])
}
#}))
print(head(bestAnnotation))

if(NROW(unique(gcSim$V1)) != NROW(bestAnnotation)) print("ERROR!") ## Sanity check....

## Get the index of each "best" annotation.
gcIndx <- match(bestAnnotation, GENCODE$V5)
tfIndx <- match(unique(gcSim$V1), TRANSCRIPT$V4)

mapped <- cbind(TRANSCRIPT[tfIndx,], GENCODE[gcIndx,])
unmapped <- TRANSCRIPT[!(TRANSCRIPT$V4 %in% unique(gcSim$V1)),]

write.table(mapped, outfileMapped, sep="\t", quote=F, row.names=F, col.names=F)
write.table(unmapped, outfileUnmapped, sep="\t", quote=F, row.names=F, col.names=F)


