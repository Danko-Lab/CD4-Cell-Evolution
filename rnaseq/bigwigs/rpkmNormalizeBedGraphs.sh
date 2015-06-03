#!/usr/bin/bash
CHINFO=../../hg19.chromInfo

##
## Normalize to RPKM for visulization.
function getCountsBw {
        echo $1
        bigWigToBedGraph $1.plus.BigWig.bw tmp.plus.bedGraph
        bigWigToBedGraph $1.minus.BigWig.bw tmp.minus.bedGraph
        cat tmp.plus.bedGraph  | awk '{print $4}' | gzip >  $1.counts.gz
        cat tmp.minus.bedGraph | awk '{print $4}' | gzip >> $1.counts.gz
	R --quiet --no-save -e 'sum(as.numeric(abs(read.table("'$1'.counts.gz")$V1)))'
        rm tmp.plus.bedGraph tmp.minus.bedGraph $1.counts.gz
}

## Get post-liftOver counts.
for i in Human.C Human.D Chimp.G Chimp.I Rhesus.J Rhesus.K
do
        echo $i
        getCountsBw $i
done


## Then ... RPKM normlize for the browser.
HC=1708934173
HD=1070415686
CG=2370446583
CI=2159259543
RJ=2641355709
RK=2130683054

function getRPKM {
	echo $1
        ## Convert to bedGraph.  Swap strands.
	bigWigToBedGraph $1.plus.BigWig.bw $1.rnaseq.minus.bg ## NOTE THE STRAND SWAP!!
        bigWigToBedGraph $1.minus.BigWig.bw $1.rnaseq.plus.bg ## NOTE STRAND SWAP.
	cat $1.rnaseq.plus.bg | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000/'"$2"'*1000000}' > $1.rnaseq.rpkm.plus.bedGraph
        cat $1.rnaseq.minus.bg | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000/'"$2"'*1000000*-1}' > $1.rnaseq.rpkm.minus.bedGraph ## Reverse here.
	bedGraphToBigWig $1.rnaseq.rpkm.plus.bedGraph $3 $1.rnaseq.rpkm.plus.bigWig
        bedGraphToBigWig $1.rnaseq.rpkm.minus.bedGraph $3 $1.rnaseq.rpkm.minus.bigWig
	rm *.bg *.bedGraph
}

getRPKM Human.C $HC ../../hg19.chromInfo
getRPKM Human.D $HD ../../hg19.chromInfo
getRPKM Chimp.G $CG ../../panTro4.chromInfo
getRPKM Chimp.I $CI ../../panTro4.chromInfo
getRPKM Rhesus.J $RJ ../../rheMac3.chromInfo
getRPKM Rhesus.K $RK ../../rheMac3.chromInfo



