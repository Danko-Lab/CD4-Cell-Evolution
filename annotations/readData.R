## 
## Read non-human primate count data, for dendrogram and expression analysis.
rowMax <- function(x) { sapply(1:NROW(x), function(i) {return(max(x[i,], na.rm=TRUE))}) }
rowMin <- function(x) { sapply(1:NROW(x), function(i) {return(min(x[i,], na.rm=TRUE))}) }

ca <- read.table("countall.tsv")
ca <- cbind(ca[,1:9], "gc18", ca[,10:NCOL(ca)])
gap <- read.table("genes.inGap")[!is.na(ca[,11]),] ## Read gap data.
ca <- ca[!is.na(ca[,11]),]

## Get pause counts.
ps <- read.table("countpause.tsv")
#ps[,7] <- paste(ps[,7],"_PauseSite", sep="")
ps <- cbind(ps[,1:6], "PauseSite", paste(ps[,4],"_PauseSite",sep=""), ps[,7], "ps", ps[,8:NCOL(ps)])
colnames(ps) <- colnames(ca)

## Get dREG counts.
ts <- read.table("counttss.tsv")
ts[,5] <- rowMax(ts[7:12])

## Rough classes...
stab <- rowMax(ts[,17:18])
dist <- ts[,13]

class <- rep("tss", NROW(ts)) ## tss is then unclassified as a promoter or enhancer
class[stab < 0.1 & dist < 500]  <- "Prox_Stab" ## Clearly protein coding promoter
class[stab > 0.1  & dist > 10000] <- "Dist_UnSt" ## Clearly distal enhancer
summary(as.factor(class))

ts <- cbind(ts[,1:6], class, ts[,c(4,20,19,21:NCOL(ts))])
#ts <- cbind(ts[,1:9], "tss", ts[,10:NCOL(ts)])
colnames(ts) <- colnames(ca)

## Join em
#ca <- rbind(ca, ps, ts)
ca <- rbind(ca, ts)  ## Don't use pause sites here!! Analyze these separately.

## Rename.
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", "annot_type",
                        "Jurkat", "Human 1", "Human 2", "Human 4", "Chimp 3", "Chimp 4", "Chimp 5", "R. Macaque 2", "R. Macaque 3", "R. Macaque 4",
                        "PI Jurkat ", "PI Human 1", "PI Human 2", "PI Human 4", "PI Chimp 3", "PI Chimp 4", "PI Chimp 5", "PI R. Macaque 2", "PI R. Macaque 3", "PI R. Macaque 4",
                        "K562", "GM12878", "IMR90", "Mouse", "Rat")
Condition <- as.factor(c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
                                                "PI", "PI", "PI", "PI", "PI", "PI", "PI", "PI", "PI", "PI", "U", "U", "U", "U", "U"))
Species  <- as.factor(c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp", "RMacaque", "RMacaque", "RMacaque",
                                                "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp", "RMacaque", "RMacaque", "RMacaque", "Human", "Human", "Human", "Mouse", "Rat"))
Labels <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                                                "Jurkat u", "Human u", "Human u", "Human u", "Chimp u", "Chimp u", "Chimp u",  "Rhe Mac u", "Rhe Mac u", "Rhe Mac u",
                                                expression(paste("Jurkat ", pi, sep="")),
                                                expression(paste("Human ", pi, sep="")),
                                                expression(paste("Human ", pi, sep="")),
                                                expression(paste("Human ", pi, sep="")),
                                                expression(paste("Chimp ", pi, sep="")),
                                                expression(paste("Chimp ", pi, sep="")),
                                                expression(paste("Chimp ", pi, sep="")),
                                                expression(paste("R Mac ", pi, sep="")),
                                                expression(paste("R Mac ", pi, sep="")),
                                                expression(paste("R Mac ", pi, sep="")),
                                                "K562", "B-cell", "IMR90" , "Mouse u", "Rat u")

## Cleanup...
ca <- ca[!is.na(ca[,11]),] ## Removes those which are not mappable/ orthologues in at least one species.
ca <- ca[grep("random", ca$chrom, invert=TRUE),]


## Used for getting useful subsets of the data.
indx.all <- c(11:35)## ALL
indx.unt <- c(11:20,31:35)## ONLY UNTREATED.
indx.good <- c(11:35) ## "Good?!"  Remove M1-PI

# Get RPKM
rpkm_df <- as.matrix(ca[,indx.good]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI
for(i in 1:NCOL(rpkm_df)) rpkm_df[,i] <- 1000*(rpkm_df[,i]+0)/sum(rpkm_df[,i]) *1000/(ca[,"mapSize"])

