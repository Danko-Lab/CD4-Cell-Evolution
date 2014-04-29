args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

## Read in main E2reg file.
E2 <- read.table(infile,header=F)
TYPE <- rep("NA", NROW(E2))
CLASS <- rep("NA", NROW(E2))

## Read in annotation files...
PATH <- "/usr/projects/GROseq/NHP/tu_caller/AnnotatePredictions/GENCODE/"

## GENCODE Annotations
refGenes <- read.table(paste(PATH,"GENCODE.Overlap.bed", sep=""))
iRG <- match(E2[[4]], refGenes[grep("protein_coding|_gene",refGenes$V15),4])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "PROTEIN_CODING"
CLASS[!is.na(iRG)] <- as.character(refGenes[grep("protein_coding|_gene",refGenes$V15),][iRG[!is.na(iRG)],12])
sum(TYPE == "PROTEIN_CODING")

iRG <- match(E2[[4]], refGenes[grep("RNA|processed_transcript|sense_overlapping|sense_intronic|3prime_overlapping_ncrna",refGenes$V15),4])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "RNA"
CLASS[!is.na(iRG)] <- as.character(refGenes[grep("RNA|processed_transcript|sense_overlapping|sense_intronic|3prime_overlapping_ncrna",refGenes$V15),][iRG[!is.na(iRG)],12])
sum(TYPE == "RNA")

iRG <- match(E2[[4]], refGenes[grep("antisense",refGenes$V15),4])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "ANTISENSE"
CLASS[!is.na(iRG)] <- as.character(refGenes[grep("antisense",refGenes$V15),][iRG[!is.na(iRG)],12])
sum(TYPE == "ANTISENSE")

iRG <- match(E2[[4]], refGenes[grep("pseudogene|GERST_PG",refGenes$V15),4])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "PSEUDOGENE+REP"
CLASS[!is.na(iRG)] <- as.character(refGenes[grep("pseudogene|GERST_PG",refGenes$V15),][iRG[!is.na(iRG)],12])
sum(TYPE == "PSEUDOGENE+REP")

## RNA Genes
rnaGenes <- read.table(paste(PATH,"rnaGene.Overlap.bed", sep=""))
iRG <- match(E2[[4]], rnaGenes[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "RNA"
CLASS[!is.na(iRG)] <- as.character(rnaGenes[iRG[!is.na(iRG)],10])
sum(TYPE == "RNA")

## Divergent transcripts
divTrans <- read.table(paste(PATH,"divergent.Overlap.bed", sep=""))
iRG <- match(E2[[4]], divTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "UPS_ANTISENSE" ## Upstream antisense (seems to be the prefered term).
CLASS[!is.na(iRG)] <- as.character(divTrans[iRG[!is.na(iRG)],10])
sum(TYPE == "UPS_ANTISENSE")

## Antisense Transc.
ansTrans <- read.table(paste(PATH,"as-refGene.Overlap.bed", sep=""))
iRG <- match(E2[[4]], ansTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "ANTISENSE"
CLASS[!is.na(iRG)] <- as.character(ansTrans[iRG[!is.na(iRG)],11])
sum(TYPE == "ANTISENSE")

## Repeat Transcription.
repTrans <- read.table(paste(PATH,"repeat.Overlap.bed", sep=""))
iRG <- match(E2[[4]], repTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "PSEUDOGENE+REP"
CLASS[!is.na(iRG)] <- as.character(repTrans[iRG[!is.na(iRG)],10])
sum(TYPE == "PSEUDOGENE+REP")

## Sense refgene bad match.
senseGeneBad <- read.table(paste(PATH, "refGeneBadMatch.Overlap.bed", sep=""))
iRG <- match(E2[[4]], senseGeneBad[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "GENE_BadMatch"
CLASS[!is.na(iRG)] <- as.character(senseGeneBad[iRG[!is.na(iRG)],11])
sum(TYPE == "GENE_BadMatch")

## Antisense refgene bad match.
asGeneBad    <- read.table(paste(PATH, "refGeneASBadMatch.Overlap.bed", sep=""))
iRG <- match(E2[[4]], asGeneBad[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "AS_BadMatch"
CLASS[!is.na(iRG)] <- as.character(asGeneBad[iRG[!is.na(iRG)],11])
sum(TYPE == "AS_BadMatch")

## Intergenic transcription
othTrans <- read.table(paste(PATH,"refASBadMatch.NoOverlap.bed", sep=""))
iRG <- match(E2[[4]], othTrans[[4]])
if(sum(TYPE[!is.na(iRG)] == "NA") != sum(!is.na(iRG))) print("POSSIBLE ERROR!") ## ERROR CHECK!
TYPE[!is.na(iRG)] <- "INTERGENIC"
CLASS[!is.na(iRG)] <- as.character(othTrans[iRG[!is.na(iRG)],4])
sum(TYPE == "INTERGENIC")

AN <- cbind(E2, TYPE= as.character(TYPE), CLASS= as.character(CLASS))
write.table(AN, outfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

###########
# Genes coated with rna...
##
## For these, I will have to make new columns.
Eout <- NULL
TYPE <- NULL
CLASS <- NULL

rnaGenes <- read.table(paste(PATH,"rnaGene.Overlap.bed", sep=""))
iRG <- match(rnaGenes[[4]], E2[[4]])

## Ensure that there are no duplicates...
rID <- paste(rnaGenes[[7]],rnaGenes[[8]],rnaGenes[[12]],sep="")
if(NROW(rnaGenes) != NROW(unique(rID))) {
  print("WARNING: POTENTIAL DUPLICATE RNAs IN DATASET: Distinct RNA Genes.")
}

for(i in c(1:NROW(iRG))) {
    Eout <- rbind(Eout, E2[iRG[i],])
    TYPE <- c(TYPE, "RNA")
    CLASS <- c(CLASS, as.character(rnaGenes[i,10]))
}

ArnaGenes <- read.table(paste(PATH,"rnaGene.Genic.Overlap.bed", sep=""))
iRG <- match(ArnaGenes[[4]], E2[[4]])

## Ensure that there are no duplicates...
rID <- paste(ArnaGenes[[7]],ArnaGenes[[8]],ArnaGenes[[12]],sep="")
if(NROW(ArnaGenes) != NROW(unique(rID))) {
  print("WARNING: POTENTIAL DUPLICATE RNAs IN DATASET: Genic Overlap.")
}

for(i in c(1:NROW(iRG))) {
    Eout <- rbind(Eout, E2[iRG[i],])
    TYPE <- c(TYPE, "RNA_withGene")
    CLASS <- c(CLASS, as.character(ArnaGenes[i,10]))
}

AN <- cbind(Eout, TYPE= as.character(TYPE), CLASS= as.character(CLASS))
write.table(AN, paste(outfile,".RNA", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

