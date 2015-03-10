#!/usr/bin/bash

## Combine data using bedtools unionbg.

## Human -U
bigWigToBedGraph Human.C.plus.BigWig.bw  H-C.plus.bedGraph
bigWigToBedGraph Human.C.minus.BigWig.bw H-C.minus.bedGraph
bigWigToBedGraph Human.D.plus.BigWig.bw  H-D.plus.bedGraph
bigWigToBedGraph Human.D.minus.BigWig.bw H-D.minus.bedGraph

bedtools unionbedg -i H-C.plus.bedGraph H-D.plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*($4+$5)}' > H.plus.bedGraph
bedtools unionbedg -i H-C.minus.bedGraph H-D.minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4+$5}' > H.minus.bedGraph

CHINFO=../../hg19.chromInfo
bedGraphToBigWig H.plus.bedGraph $CHINFO H.plus.bw
bedGraphToBigWig H.minus.bedGraph $CHINFO H.minus.bw

## Chimp -U
#bigWigToBedGraph Chimp.F.plus.BigWig.bw  C-F.plus.bedGraph
#bigWigToBedGraph Chimp.F.minus.BigWig.bw C-F.minus.bedGraph
bigWigToBedGraph Chimp.G.plus.BigWig.bw  C-G.plus.bedGraph
bigWigToBedGraph Chimp.G.minus.BigWig.bw C-G.minus.bedGraph
#bigWigToBedGraph Chimp.H.plus.BigWig.bw  C-H.plus.bedGraph
#bigWigToBedGraph Chimp.H.minus.BigWig.bw C-H.minus.bedGraph
bigWigToBedGraph Chimp.I.plus.BigWig.bw  C-I.plus.bedGraph
bigWigToBedGraph Chimp.I.minus.BigWig.bw C-I.minus.bedGraph

bedtools unionbedg -i C-G.plus.bedGraph C-I.plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*($4+$5)}' > C.plus.bedGraph
bedtools unionbedg -i C-G.minus.bedGraph C-I.minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4+$5}' > C.minus.bedGraph

CHINFO=../../panTro4.chromInfo
bedGraphToBigWig C.plus.bedGraph $CHINFO C.plus.bw
bedGraphToBigWig C.minus.bedGraph $CHINFO C.minus.bw

## Rhesus -U
bigWigToBedGraph Rhesus.J.plus.BigWig.bw  R-J.plus.bedGraph
bigWigToBedGraph Rhesus.J.minus.BigWig.bw R-J.minus.bedGraph
bigWigToBedGraph Rhesus.K.plus.BigWig.bw  R-K.plus.bedGraph
bigWigToBedGraph Rhesus.K.minus.BigWig.bw R-K.minus.bedGraph

bedtools unionbedg -i R-J.plus.bedGraph R-K.plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*($4+$5)}' > R.plus.bedGraph
bedtools unionbedg -i R-J.minus.bedGraph R-K.minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-$4+$5}' > R.minus.bedGraph

CHINFO=../../rheMac3.chromInfo
bedGraphToBigWig R.plus.bedGraph $CHINFO R.plus.bw
bedGraphToBigWig R.minus.bedGraph $CHINFO R.minus.bw

rm *.bedGraph

##
## Swap strands.
mv H.plus.bw H-U.rnaseq.minus.bw
mv H.minus.bw H-U.rnaseq.plus.bw

mv C.plus.bw C-U.rnaseq.minus.bw
mv C.minus.bw C-U.rnaseq.plus.bw

mv R.plus.bw M-U.rnaseq.minus.bw
mv R.minus.bw M-U.rnaseq.plus.bw

##
## Normalize to RPKM for visulization.
function getCountsBw {
        echo $1
        bigWigToBedGraph $1.rnaseq.plus.bw tmp.plus.bedGraph
        bigWigToBedGraph $1.rnaseq.minus.bw tmp.minus.bedGraph
        cat tmp.plus.bedGraph  | awk '{print $4}' | gzip >  $1.counts.gz
        cat tmp.minus.bedGraph | awk '{print $4}' | gzip >> $1.counts.gz
	R --quiet --no-save -e 'sum(as.numeric(abs(read.table("'$1'.counts.gz")$V1)))'
        rm tmp.plus.bedGraph tmp.minus.bedGraph $1.counts.gz
}

## Get post-liftOver counts.
for i in H-U C-U M-U
do
        echo $i
        getCountsBw $i
done


## Then ... RPKM normlize for the browser.
HCOUNTS=3318960093
CCOUNTS=5006443463
MCOUNTS=2958335679


function getRPKM {
	echo $1
	bigWigToBedGraph $1.rnaseq.plus.bw $1.rnaseq.plus.bg
        bigWigToBedGraph $1.rnaseq.minus.bw $1.rnaseq.minus.bg
	cat $1.rnaseq.plus.bg | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000/'"$2"'*1000000}' > $1.rnaseq.rpkm.plus.bedGraph
        cat $1.rnaseq.minus.bg | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000/'"$2"'*1000000}' > $1.rnaseq.rpkm.minus.bedGraph
	bedGraphToBigWig $1.rnaseq.rpkm.plus.bedGraph $3 $1.rnaseq.rpkm.plus.bigWig
        bedGraphToBigWig $1.rnaseq.rpkm.minus.bedGraph $3 $1.rnaseq.rpkm.minus.bigWig
	rm *.bg *.bedGraph
}

getRPKM H-U $HCOUNTS ../../../hg19.chromInfo
getRPKM C-U $HCOUNTS ../../../panTro4.chromInfo
getRPKM M-U $HCOUNTS ../../../rheMac3.chromInfo


