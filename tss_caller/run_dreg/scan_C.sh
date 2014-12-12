#$ -S /bin/bash
#$ -cwd
#$ -N dREG.chimp
#$ -o dREG.chimp.out.$JOB_ID
#$ -j y
#$ -pe bscb 16
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=48:00:00
#$ -q long_term.q

PREFIX=C
DATADIR=AllData/All_Merge
TH=0.8

## Copy files to scratch space (/workdir and /SSD).
STARTDIR=`pwd`
SCRATCH=/SSD/cgd24_tssDetector_$PREFIX
mkdir $SCRATCH
cp /bscb/bscb07/cgd24/projects/NHP/tss_caller/run_dreg/scan_nhp.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/train_svm/asvm.intersDNase.getTrainSet.RData $SCRATCH
#cp /home/cgd24/projects/tss_detector/train_svm_nhp/cd4.dnase1.adbn.RData $SCRATCH
#cp /home/cgd24/projects/tss_detector/train_svm/asvm.RData $SCRATCH ## 
cp /bscb/bscb07/cgd24/projects/NHP/$DATADIR/$PREFIX*us.bw $SCRATCH
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
