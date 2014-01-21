## 
## Read non-human primate count data, for dendrogram and expression analysis.

ca <- read.table("countall.tsv")
ca <- cbind(ca[,1:9], "gc18", ca[,10:NCOL(ca)])

## Get pause counts.
ps <- read.table("countpause.tsv")
ps[,7] <- paste(ps[,7],"_PauseSite", sep="")
ps <- cbind(ps[,1:9], "ps", ps[,10:NCOL(ps)])
colnames(ps) <- colnames(ca)

## Get dREG counts.
ts <- read.table("counttss.tsv")
ts <- cbind(ts[,1:6], "PromEnh",  ts[,c(4,13:36)])
ts <- cbind(ts[,1:9], "tss", ts[,10:NCOL(ts)])
colnames(ts) <- colnames(ca)

## Join em
ca <- rbind(ca, ps, ts)

## Rename.
names(ca) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "type", "mgi", "mapSize", "annot_type",
                                                "Jurkat", "Human 1", "Human 2", "Human 3", "Chimp 2", "Chimp 3", "Chimp 4", "R. Macaque 1", "R. Macaque 2", "R. Macaque 3",
                                                "PI Jurkat ", "PI Human 1", "PI Human 2", "PI Human 3", "PI Chimp 2", "PI Chimp 3", "PI Chimp 4", "PI R. Macaque 1", "PI R. Macaque 2", "PI R. Macaque 3",
                                                "K562", "GM12878", "IMR90")
Condition <- as.factor(c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                "U", "U", "U", "U", "U", "U", "U", "U", "U", "U",
                                                "PI", "PI", "PI", "PI", "PI", "PI", "PI", "PI", "PI", "PI", "U", "U", "U"))
Species  <- as.factor(c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", 
                                                "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp", "RMacaque", "RMacaque", "RMacaque",
                                                "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp", "RMacaque", "RMacaque", "RMacaque", "Human", "Human", "Human"))
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
                                                "K562", "B-cell", "IMR90")

ca <- ca[!is.na(ca[,11]),]


## Used for getting useful subsets of the data.
indx.all <- c(11:33)## ALL
indx.unt <- c(11:14,16:20,31:33)## ONLY UNTREATED, GOOD Remove C2-U
#indx <- c(11:12,14:22,24,26) ## "Good?!"  Remove H2-U, C2-PI, M1-PI
indx.good <- c(11:14,16:24,26:27,29:33) ## "Good?!"  Remove C2-U+PI, M1-PI


