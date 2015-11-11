library(parallel);
options("scipen"=100, "digits"=4)

load("union-test.rdata");

get_chromosome_size <- function(file.twoBit)
{
	file.tmp <- tempfile();

	err_code <- system(paste("twoBitInfo", file.twoBit, file.tmp, sep=" "));
	chromInfo <- read.table( file.tmp );
	unlink( file.tmp );

	if( err_code!=0 ) cat("! Failed to call the twoBitInfo command\n");

	return(chromInfo);
}
get_complement<-function( bed, file.twoBit )
{
	tmp.bed <- tempfile();
	write.table(bed, file=tmp.bed, quote=F, sep="\t", row.names=FALSE, col.names=FALSE,);

	tmp.genome <- tempfile();
	write.table( get_chromosome_size(file.twoBit), file=tmp.genome, quote=F, sep="\t", row.names=FALSE, col.names=FALSE,);

	pipe.cmd <- paste("bedtools complement -i", tmp.bed, " -g ", tmp.genome );
	df.comp <- read.table( pipe(pipe.cmd), header = F );

	return(df.comp);
}

get_intersect <- function( bedA, bedB )
{
	tmp.bedA <- tempfile();
	write.table(bedA, file=tmp.bedA, quote=F, sep="\t", row.names=FALSE, col.names=FALSE);

	tmp.bedB <- tempfile();
	write.table(bedB, file=tmp.bedB, quote=F, sep="\t", row.names=FALSE, col.names=FALSE);

	pipe.cmd<- paste( "bedtools intersect -a ", tmp.bedA, "-b", tmp.bedB );
	df.bed <- try( read.table( pipe(pipe.cmd), header = F ) );
	
	if(class(df.bed)=="try-error")
		return(c());

	unlink( tmp.bedA );
	unlink( tmp.bedB );
	
	return(df.bed);
}


file.twoBit.hg19 <- "/local/storage/data/hg19/hg19.2bit";

##Differences in DNA sequences (bed file C):
fileC <- "/local/storage/projects/NHP/dna_sequence/hg19.diff.bed.starch"

## Differences in dREG (bed file B):
fileB <- "/local/storage/projects/NHP/enhancer/hg19.gain.loss.bed"

T.B <- read.table(fileB, sep="\t", header=F);
T.C <- read.table(pipe(paste("unstarch ", fileC)), sep="\t", header=F);

T.notB <- get_complement( T.B,  file.twoBit.hg19 );
T.notC <- get_complement( T.C,  file.twoBit.hg19 );

T.BC   <- get_intersect( T.B,    T.C)
T.nBnC <- get_intersect( T.notB, T.notC)
T.nBC  <- get_intersect( T.notB, T.C)
T.BnC  <- get_intersect( T.B,    T.notC)

get_fisher_test <- function( bed )
{
	t11 <- NROW( get_intersect( bed,  T.BC) )
	t12 <- NROW( get_intersect( bed,  T.nBC)  )
	t21 <- NROW( get_intersect( bed,  T.BnC)  )
	t22 <- NROW( get_intersect( bed,  T.nBnC)  )

	pv <- fisher.test( rbind(c(t11,t12), c(t21,t22)) )
	return(c(pv$p.value, t11, t12, t21, t22));
}

#u.list <- mclapply( 1:length(tf.human.u$result), function(i){
u.list <- mclapply( 1:200, function(i){
cat(i, "\n");
	p <- try(get_fisher_test( tf.human.u$result[[i]] ));
	p;
}, mc.cores=21);

u.summary <- do.call(rbind, u.list)


pi.list <- mclapply( 1:length(tf.human.pi$result), function(i){
	get_fisher_test( tf.human.pi$result[[i]] );
}, mc.cores=7);

pi.summary <- do.call(rbind, pi.list)
