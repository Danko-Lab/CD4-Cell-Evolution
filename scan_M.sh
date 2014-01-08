#$ -S /bin/bash
#$ -cwd
#$ -N dREG.macaque
#$ -o dREG.macaque.out.$JOB_ID
#$ -j y
#$ -pe bscb 32
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

PREFIX=M
DATADIR=All_Merge
TH=0.83

## Copy files to scratch space (/workdir and /SSD).
STARTDIR=`pwd`
SCRATCH=/SSD/cgd24_tssDetector_$PREFIX
mkdir $SCRATCH
cp $STARTDIR/scan_nhp.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp ~/nextgen/projects/GROseq/NHP/$DATADIR/$PREFIX*.bw $SCRATCH
cd $SCRATCH

## Run R.
R --no-save --args $PREFIX < scan_nhp.R

## Merge into TSS. 
for i in `ls *.TSS.bedGraph.gz`
do
  j=`echo $i | cut -d \. -f 1,2`
  echo $j

  ## Merge into peaks.
  zcat $i | awk 'BEGIN{OFS="\t"} ($4 > '"$TH"') {print $1,$2-50,$3+51}' | sort-bed - | bedops --merge - > $j.no_score.bed

  ## Get max score in each peak. 
  zcat $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4}' | sort-bed - | bedmap --echo --max $j.no_score.bed - | sed "s/|/\t/g" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"n",$4}' | gzip > $j.bed.gz

  ## Clean temporary files.  
  rm $j.no_score.bed
done

## Copy data files back.
cp *.TSS.bedGraph.gz $STARTDIR
cp *.bed.gz $STARTDIR

## Cleanup
rm -Rf $SCRATCH
