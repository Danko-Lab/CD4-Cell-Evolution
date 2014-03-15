############################################################
##
## General summaries of changes in pausing ... compared to gene expression.

source("../lib/densScatterplot.R")
#example: densScatterplot(a$V7, a$V8, xlab="Human Unt.", ylab="Chimpanzee Unt.", main="dREG Scores")

load("../annotations/fdr.RData")

## Sanity check ... use Amean?!
cor.test(hs$Amean, log(ca$"Human 1"+1))
densScatterplot(hs$Amean, log(ca$"Human 1"+1), xlab="Limma", ylab="Raw counts")

## Compare Human and Chimp.
cor.test(ca$"Human 1", ca$"Chimp 4")
densScatterplot(log(ca$"Human 1"+1), log(ca$"Chimp 4"+1), xlab="Human 1", ylab="Chimpanzee 4", main="All TU")

## Just protein-coding genes.
cmpPauseBody <- function(sample1, sample2, epsilon=1e-4) {
 par(mfrow=c(1,2))

 print("Protein Coding")
 indx <- ca$type == "protein_coding" & ca[,sample1]>0 & ca[,sample2]>0 ## Requiring >0 reads in both samples prevents falling into an alignment gap.
 print(cor.test(ca[indx,sample1]/ca$mapSize[indx], ca[indx,sample2]/ca$mapSize[indx], method="spearman"))
 densScatterplot(log(ca[indx,sample1]/ca$mapSize[indx]+epsilon), 
		log(ca[indx,sample2]/ca$mapSize[indx]+epsilon), 
		xlab=sample1, ylab=sample2, main="Protein-coding genes")
 abline(0,1)

 print("Pause Sites")
 indx <- ca$type == "protein_coding_PauseSite" & ca[,sample1]>0 & ca[,sample2]>0
 print(cor.test(ca[indx,sample1]/ca$mapSize[indx], ca[indx,sample2]/ca$mapSize[indx], method="spearman"))
 densScatterplot(log(ca[indx,sample1]/ca$mapSize[indx]+epsilon), 
		log(ca[indx,sample2]/ca$mapSize[indx]+epsilon), 
		xlab=sample1, ylab=sample2, main="Pause sites")
 abline(0,1)
}

cmpPauseBody("Human 1", "Human 2")
cmpPauseBody("Human 1", "Human 3")

cmpPauseBody("Human 1", "Chimp 4")
cmpPauseBody("Human 1", "R. Macaque 3")

cmpPauseBody("Human 1", "R. Macaque 2")
cmpPauseBody("Human 1", "R. Macaque 1")

## Compute some sort of RMSD for changes in pausing levels.
## Be sure to quantile normalize.

############################################################
##
## Get pausing indices.

indxPause <- ca$type=="protein_coding_PauseSite"
indxBody  <- ca$type=="protein_coding"

body <- ca[indxBody,][ca[indxBody,"name"] %in% ca[indxPause,"name"],]
pause <- pause <- ca[indxPause,][match(as.character(body[,"name"]), as.character(ca[indxPause,"name"])),] 
stopifnot(sum(body$name == pause$name) == NROW(pause)) ## SANTIY CHECK.

indx <- pause[,"Human 1"]>0 & pause[,"Chimp 4"]>0 & pause[,"R. Macaque 3"]>0 ## Unnecesarily conservative.
PI <- (((pause[,c(11:NCOL(pause))]+1)/pause[,"mapSize"]) / ((body[,c(11:NCOL(body))]+1)/body[,"mapSize"]))[indx,]

cor.test(PI[,"Human 1"], PI[,"Human 2"], method="spearman")
cor.test(PI[,"Human 1"], PI[,"Chimp 4"], method="spearman")
cor.test(PI[,"Human 1"], PI[,"R. Macaque 3"], method="spearman")

## Plot...
par(mfrow=c(1,3))
densScatterplot(log(PI[,"Human 1"]), log(PI[,"Human 2"]), xlab="Human 1", ylab="Human 2")
densScatterplot(log(PI[,"Human 1"]), log(PI[,"Chimp 4"]), xlab="Human 1", ylab="Chimp 4")
densScatterplot(log(PI[,"Human 1"]), log(PI[,"R. Macaque 3"]), xlab="Human 1", ylab="R. Macaque 3")

############################################################



