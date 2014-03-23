#$ -S /bin/bash
#$ -cwd
#$ -N train_dbn
#$ -o train_dbn.out.$JOB_ID
#$ -j y
#$ -pe bscb 8
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

STARTDIR=`pwd`

## Copy files to scratch space (/workdir and /SSD).
SCRATCH=/SSD/cgd24_traindbn_cd4/
mkdir $SCRATCH
cd $SCRATCH

cp ~/nextgen/projects/GROseq/NHP/tss_caller/retrain_dreg/train_dbn.R $SCRATCH ## 
cp /home/cgd24/projects/tss_detector/data/GencodeMerge.IntersectOpStrand.bed $SCRATCH ## Gene overlap files

## CD4
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/cd4/dnase1fp/dnase1.peaks_peaks.narrowPeak $SCRATCH ## DNAse-1
cp ~/nextgen/projects/GROseq/NHP/AllData/All_Merge/*.bw $SCRATCH ## PRO-seq
cp /home/cgd24/nextgen/data/GROseq.parser/hg19/cd4/chromhmm/CD4.chromHMM.Ernst2010.hg19.Prom.Enh.bed $SCRATCH ## Chromatin marks.
## Run R.
R --no-save < train_dbn.R

## Copy data files back.
cp *.bedGraph $STARTDIR
cp *.RData $STARTDIR
cp *.pdf $STARTDIR

## Cleanup
rm -Rf $SCRATCH

