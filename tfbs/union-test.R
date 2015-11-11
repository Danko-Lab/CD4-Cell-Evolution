library(rtfbsdb);
library(parallel);
library(data.table);
options("scipen"=100, "digits"=4);
ncores=42;

	
merge_bed_u_hg19<-function()
{
	bed.dreg.chimp.u.hg19   <- crossmap_bed_to_species( file.dreg.chimp.u,   file.chain.chimp.human);
	bed.dreg.rhesus.u.hg19  <- crossmap_bed_to_species( file.dreg.rhesus.u,  file.chain.rhesus.human);
	bed.dreg.human.u  <- read.table(file.dreg.human.u, header=F);
	bed.dreg.u  <- rbind( bed.dreg.human.u, bed.dreg.chimp.u.hg19, bed.dreg.rhesus.u.hg19);

	n.wrong.u <- length(which(bed.dreg.u[,2]>=bed.dreg.u[,3]))
	if(n.wrong.u>0)
	{
		cat("START>END", length(which(bed.dreg.u[,2]>=bed.dreg.u[,3])), "\n");
		bed.dreg.u <- bed.dreg.u [ -which(bed.dreg.u[,2]>=bed.dreg.u[,3]), ]
	}

	bed.dreg.u  <- bedtools_merge_bed(bed.dreg.u);

	return(bed.dreg.u);
}

merge_bed_pi_hg19<-function()
{
	bed.dreg.chimp.pi.hg19  <- crossmap_bed_to_species( file.dreg.chimp.pi,  file.chain.chimp.human);
	bed.dreg.rhesus.pi.hg19 <- crossmap_bed_to_species( file.dreg.rhesus.pi, file.chain.rhesus.human);
	bed.dreg.human.pi <- read.table(file.dreg.human.pi, header=F);
	bed.dreg.pi <- rbind( bed.dreg.human.pi, bed.dreg.chimp.pi.hg19, bed.dreg.rhesus.pi.hg19);

	n.wrong.pi <- length(which(bed.dreg.pi[,2]>=bed.dreg.pi[,3]))
	if(n.wrong.pi>0)
	{
		cat("START>END", length(which(bed.dreg.pi[,2]>=bed.dreg.pi[,3])), "\n");
		bed.dreg.pi <- bed.dreg.pi [ -which(bed.dreg.pi[,2]>=bed.dreg.pi[,3]), ]
	}
	bed.dreg.pi <- bedtools_merge_bed(bed.dreg.pi);
	
	return(bed.dreg.pi);
}

bedtools_merge_bed<-function(bed)
{
	tmp.bed <- tempfile();
	write.table(bed, file=tmp.bed, quote=F, sep="\t", row.names=FALSE, col.names=FALSE,);

	pipe.cmd <- paste("sort-bed ", tmp.bed, " | bedtools merge  -i - ");
	df.merge <- read.table( pipe(pipe.cmd), header = F );

	return(df.merge);
}

crossmap_bed_to_species <- function( bed, file.chain)
{
	file.bed <- bed ;
	if (class(bed)=="data.frame")
	{
		file.bed <- tempfile();
		write.table(bed, file=file.bed, sep="\t", quote=F, col.names=F, row.names=F);
	}
	
	file.tmp <- tempfile();
	system( paste("CrossMap.py bed ", file.chain, file.bed, file.tmp) )
	df <- read.table( file.tmp, header=F, sep="\t");	
	
	return(df);
}

crossmap_tfbs_result <- function(tfbs.finding, file.chain)
{
	bed.all <- as.data.frame( rbindlist(tfbs.finding$result) );
	bed.change <- crossmap_bed_to_species( bed.all, file.chain);
	
	for(i in 1:length(tfbs.finding$result))
		tfbs.finding$result[[i]] <- bed.change[ which( bed.change[,4] == as.character(tfbs.finding$result[[i]][1,4]) ), ]	
	
	return(tfbs.finding);
}


# U = A +B + C
get_union_bed3<-function( bed1, bed2, bed3 )
{
	bed1.str <- bed2.str <- bed3.str <- c();
	
	if (!is.null(bed1) && NROW(bed1)>0 ) bed1.str <- paste(bed1[,1], ":", bed1[,2], ":", bed1[,3], sep="");
	if (!is.null(bed2) && NROW(bed2)>0 ) bed2.str <- paste(bed2[,1], ":", bed2[,2], ":", bed2[,3], sep="");
	if (!is.null(bed3) && NROW(bed3)>0 ) bed3.str <- paste(bed3[,1], ":", bed3[,2], ":", bed3[,3], sep="");

	bed.str <- unique(c(bed1.str, bed2.str, bed3.str));
	if(is.null(bed.str))
		return(c());
		
	bed.all <- as.data.frame(do.call(rbind, strsplit(bed.str, ":")));
	
	file.bed1 <- tempfile();	
	write.table( bed.all, file=file.bed1, quote=F, row.names=F, col.names=F, sep="\t");

	pipe.cmd <- paste("sort-bed ", file.bed1, sep=" ");
	bed3.union <- try( read.table( pipe(pipe.cmd), header = F ) );
	
	if(class(bed3.union)=="try-error")
		bed3.union <- c();
	
	unlink(file.bed1);
	
	return(bed3.union);
}

write.starchbed <- function(bed, file.starch) {
	# pipe bed into starch file
	write.table(bed, file = pipe(paste("sort-bed - | starch - >", file.starch)), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	if( !file.exists(file.starch) )    
		cat("! Failed to write starch file (", file.starch, ") using the sort-bed and starch commands.\n");
}


merge_three_scanning<- function( bed1.finding, bed2.finding, bed3.finding, motif.id, tf.name )
{
	bed3.union <- get_union_bed3( bed1.finding, bed2.finding, bed3.finding );
	if(NROW(bed3.union)<1) return(bed3.union);
	
	bed3.score <- apply( bed3.union, 1, function(x){

		bed1.match <- bed2.match <- bed3.match <- c()

		if (!is.null(bed1.finding)) bed1.match <- which( as.character(x[1])==as.character(bed1.finding[,1]) & as.numeric(x[2])==bed1.finding[,2] & as.numeric(x[3])==bed1.finding[,3])
		if (!is.null(bed2.finding)) bed2.match <- which( as.character(x[1])==as.character(bed2.finding[,1]) & as.numeric(x[2])==bed2.finding[,2] & as.numeric(x[3])==bed2.finding[,3])
		if (!is.null(bed3.finding)) bed3.match <- which( as.character(x[1])==as.character(bed3.finding[,1]) & as.numeric(x[2])==bed3.finding[,2] & as.numeric(x[3])==bed3.finding[,3])
		
		score1 <- score2 <- score3 <- NA;
		
		#if(length(bed1.match)>1 || length(bed2.match)>1 || length(bed3.match)>1 ) browser();
		
		if(length(bed1.match)>0)  score1 <- bed1.finding[bed1.match[1], 5];
		if(length(bed2.match)>0)  score2 <- bed2.finding[bed2.match[1], 5];
		if(length(bed3.match)>0)  score3 <- bed3.finding[bed3.match[1], 5];
		
		return(data.frame(score1, score2, score3));
	})
	
	#bed3.score <- do.call(rbind, bed3.score);
	bed3.score <- rbindlist( bed3.score );
	
	max.score <- apply(bed3.score, 1, function(x) max(x, na.rm=T) );
	
	bed6.union <- data.frame( bed3.union, motif.id, tf.name, max.score, bed3.score );
	
	return(bed6.union);
}

merge_tf <- function(tf.human, tf.chimp, tf.rhesus)
{
	tf.all.motif <- unique(rbind( tf.human$summary[,c(1,2)], 
			tf.chimp$summary[,c(1,2)], 
			tf.rhesus$summary[,c(1,2)] ));

	result <- mclapply( 1:NROW(tf.all.motif), function(i){

			idx1 <- idx2 <- idx3 <- c();
			idx1 <- which( as.character(tf.human$summary[,1])  == as.character( tf.all.motif[i,1]) & tf.human$summary[,2] == as.character( tf.all.motif[i,2]) )
			idx2 <- which( as.character(tf.chimp$summary[,1])  == as.character( tf.all.motif[i,1]) & tf.chimp$summary[,2] == as.character( tf.all.motif[i,2]) )
			idx3 <- which( as.character(tf.rhesus$summary[,1]) == as.character( tf.all.motif[i,1]) & tf.rhesus$summary[,2] == as.character( tf.all.motif[i,2]) );

			bed1 <- bed2 <- bed3 <- NULL;
			if(length(idx1)>0) bed1 <- tf.human$result[[idx1]];
			if(length(idx2)>0) bed2 <- tf.chimp$result[[idx2]];
			if(length(idx3)>0) bed3 <- tf.rhesus$result[[idx3]];

			r <- try(merge_three_scanning(bed1, bed2, bed3, 
				as.character( tf.all.motif[i,1]), 
				as.character( tf.all.motif[i,2]) ) );

			cat("TF=", i, "\t", as.character( tf.all.motif[i,1]), 
				"\t", as.character( tf.all.motif[i,2]), "\t", idx1, "\t", idx2, "\t", idx3);

			if( class(r)!="try-error" )		
				cat( "\tR=", NROW(r), "\n")
			else
			{
				show(r);
				r <- as.character(r);
			}
			
			return(r);
			
		}, mc.cores = ncores );

	return(list(summary=tf.all.motif, result=result));
}


## basic data
file.hg19.twoBit    <-  "/local/storage/data/hg19/hg19.2bit";
file.panTro4.twoBit <-  "/local/storage/data/2bit/panTro4.2bit";
file.rheMac3.twoBit <-  "/local/storage/data/2bit/rheMac3.2bit";
file.gencode.gtf    <-  "/local/storage/data/gencode/gencode.v19.annotation.gtf";

## dREG-HD: 
file.dreg.human.u   <- "/local/storage/projects/NHP/dREG_HD/H-U_dREG_HD.bed"
file.dreg.human.pi  <- "/local/storage/projects/NHP/dREG_HD/H-PI_dREG_HD.bed"
file.dreg.chimp.u   <- "/local/storage/projects/NHP/dREG_HD/C-U_dREG_HD.bed"
file.dreg.chimp.pi  <- "/local/storage/projects/NHP/dREG_HD/C-PI_dREG_HD.bed"
file.dreg.rhesus.u  <- "/local/storage/projects/NHP/dREG_HD/M-U_dREG_HD.bed"
file.dreg.rhesus.pi <- "/local/storage/projects/NHP/dREG_HD/M-PI_dREG_HD.bed"

## Chains (for crossmap): 
## Rhesus -> Human: 
file.chain.rhesus.human <- "/local/storage/projects/NHP/makeRecipBest/hg19.rheMac3/rheMac3.hg19.rbest.chain.gz"
## Chimp -> Human: 
file.chain.chimp.human  <- "/local/storage/projects/NHP/makeRecipBest/hg19.panTro4/panTro4.hg19.rbest.chain.gz"
## Human -> Rhesus: 
file.chain.human.rhesus <- "/local/storage/projects/NHP/makeRecipBest/hg19.rheMac3/axtChain/hg19.rheMac3.rbest.chain.gz"
## Human -> Chimp: 
file.chain.human.chimp  <- "/local/storage/projects/NHP/makeRecipBest/hg19.panTro4/axtChain/hg19.panTro4.rbest.chain.gz"

## PRO-seq data: 
## Human: 
file.human.u.minus.bw   <- "/local/storage/projects/NHP/AllData/All_Merge/H-U_minus.hg19.bw"
file.human.u.plus.bw    <- "/local/storage/projects/NHP/AllData/All_Merge/H-U_plus.hg19.bw" 
file.human.pi.minus.bw  <- "/local/storage/projects/NHP/AllData/All_Merge/H-PI_minus.hg19.bw"
file.human.pi.plus.bw   <- "/local/storage/projects/NHP/AllData/All_Merge/H-PI_plus.hg19.bw"

## Chimp: 
file.chimp.u.minus.bw   <- "/local/storage/projects/NHP/AllData/All_Merge/C-U_minus.hg19.bw"
file.chimp.u.plus.bw    <- "/local/storage/projects/NHP/AllData/All_Merge/C-U_plus.hg19.bw"
file.chimp.pi.minus.bw  <- "/local/storage/projects/NHP/AllData/All_Merge/C-PI_minus.hg19.bw"
file.chimp.pi.plus.bw   <- "/local/storage/projects/NHP/AllData/All_Merge/C-PI_plus.hg19.bw"

## Rhesus Macaque: 
file.rhesus.u.minus.bw  <- "/local/storage/projects/NHP/AllData/All_Merge/M-U_minus.hg19.bw"
file.rhesus.u.plus.bw   <- "/local/storage/projects/NHP/AllData/All_Merge/M-U_plus.hg19.bw"
file.rhesus.pi.minus.bw <- "/local/storage/projects/NHP/AllData/All_Merge/M-PI_minus.hg19.bw"
file.rhesus.pi.plus.bw  <- "/local/storage/projects/NHP/AllData/All_Merge/M-PI_plus.hg19.bw"

db <- CisBP.extdata( "human" );
tfs <- tfbs.createFromCisBP( db ); 

tfs.human.u  <- tfbs.selectExpressedMotifs( tfs, file.hg19.twoBit , 
				file.gencode.gtf, 
				file.bigwig.plus = file.human.u.plus.bw, 
				file.bigwig.minus = file.human.u.minus.bw, 
				seq.datatype ="GRO-seq", 
				pvalue.threshold = 0.05,
				ncores = ncores); 
				
tfs.human.pi <- tfbs.selectExpressedMotifs( tfs, file.hg19.twoBit , 
				file.gencode.gtf, 
				file.bigwig.plus = file.human.pi.plus.bw, 
				file.bigwig.minus = file.human.pi.minus.bw, 
				seq.datatype = "GRO-seq", 
				pvalue.threshold = 0.05,
				ncores = ncores); 

tfs.chimp.u  <- tfbs.selectExpressedMotifs( tfs, file.hg19.twoBit , 
				file.gencode.gtf, 
				file.bigwig.plus = file.chimp.u.plus.bw, 
				file.bigwig.minus = file.chimp.u.minus.bw, 
				seq.datatype ="GRO-seq", 
				pvalue.threshold = 0.05,
				ncores = ncores); 

tfs.chimp.pi <- tfbs.selectExpressedMotifs( tfs, file.hg19.twoBit , 
				file.gencode.gtf, 
				file.bigwig.plus = file.chimp.pi.plus.bw, 
				file.bigwig.minus = file.chimp.pi.minus.bw, 
				seq.datatype ="GRO-seq", 
				pvalue.threshold = 0.05,
				ncores = ncores); 

tfs.rhesus.u <- tfbs.selectExpressedMotifs( tfs, file.hg19.twoBit , 
				file.gencode.gtf, 
				file.bigwig.plus = file.rhesus.u.plus.bw, 
				file.bigwig.minus = file.rhesus.u.minus.bw, 
				seq.datatype = "GRO-seq", 
				pvalue.threshold = 0.05,
				ncores = ncores); 

tfs.rhesus.pi <- tfbs.selectExpressedMotifs( tfs, file.hg19.twoBit , 
				file.gencode.gtf, 
				file.bigwig.plus = file.rhesus.pi.plus.bw, 
				file.bigwig.minus = file.rhesus.pi.minus.bw, 
				seq.datatype = "GRO-seq", 
				pvalue.threshold = 0.05,
				ncores = ncores); 

bed.dreg.u.hg19    <- merge_bed_u_hg19();
bed.dreg.u.chimp   <- crossmap_bed_to_species( bed.dreg.u.hg19, file.chain.human.chimp);
bed.dreg.u.rhesus  <- crossmap_bed_to_species( bed.dreg.u.hg19, file.chain.human.rhesus);

bed.dreg.pi.hg19   <- merge_bed_pi_hg19();
bed.dreg.pi.chimp  <- crossmap_bed_to_species( bed.dreg.pi.hg19, file.chain.human.chimp);
bed.dreg.pi.rhesus <- crossmap_bed_to_species( bed.dreg.pi.hg19, file.chain.human.rhesus);

bed.dreg.u.hg19    <- bed.dreg.u.hg19 [grep("_|chrM|chrY|chrX", bed.dreg.u.hg19[,1], invert=TRUE),];
bed.dreg.u.chimp   <- bed.dreg.u.chimp [grep("_|chrM|chrY|chrX", bed.dreg.u.chimp[,1], invert=TRUE),];
bed.dreg.u.rhesus  <- bed.dreg.u.rhesus [grep("_|chrM|chrY|chrX", bed.dreg.u.rhesus[,1], invert=TRUE),];

tf.human.u         <- tfbs.scanTFsite( tfs.human.u, file.hg19.twoBit, bed.dreg.u.hg19, threshold = 7, ncores = ncores);
tf.chimp.u         <- tfbs.scanTFsite( tfs.chimp.u, file.panTro4.twoBit, bed.dreg.u.chimp, threshold = 7, ncores = ncores)
tf.chimp.u         <- crossmap_tfbs_result( tf.chimp.u, file.chain.chimp.human );
tf.rhesus.u        <- tfbs.scanTFsite( tfs.rhesus.u, file.rheMac3.twoBit, bed.dreg.u.rhesus, threshold = 7, ncores = ncores)
tf.rhesus.u        <- crossmap_tfbs_result( tf.rhesus.u, file.chain.rhesus.human );

bed.dreg.pi.hg19   <- bed.dreg.pi.hg19[grep("_|chrM|chrY|chrX", bed.dreg.pi.hg19[,1], invert=TRUE),];
bed.dreg.pi.chimp  <- bed.dreg.pi.chimp[grep("_|chrM|chrY|chrX", bed.dreg.pi.chimp[,1], invert=TRUE),];
bed.dreg.pi.rhesus <- bed.dreg.pi.rhesus[grep("_|chrM|chrY|chrX", bed.dreg.pi.rhesus[,1], invert=TRUE),];

tf.human.pi        <- tfbs.scanTFsite( tfs.human.pi, file.hg19.twoBit, bed.dreg.pi.hg19, threshold = 7, ncores = ncores);
tf.chimp.pi        <- tfbs.scanTFsite( tfs.chimp.pi, file.panTro4.twoBit, bed.dreg.pi.chimp, threshold = 7, ncores = ncores)
tf.chimp.pi        <- crossmap_tfbs_result( tf.chimp.pi, file.chain.chimp.human );
tf.rhesus.pi       <- tfbs.scanTFsite( tfs.rhesus.pi, file.rheMac3.twoBit, bed.dreg.pi.rhesus, threshold = 7, ncores = ncores)
tf.rhesus.pi       <- crossmap_tfbs_result( tf.rhesus.pi, file.chain.rhesus.human );

save.image("union-test.rdata");

tf.merge.u  <- merge_tf( tf.human.u,  tf.chimp.u,  tf.rhesus.u );
tf.merge.pi <- merge_tf( tf.human.pi, tf.chimp.pi, tf.rhesus.pi);

save.image("union-test.rdata");