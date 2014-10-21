#$ -S /bin/bash
#$ -cwd
#$ -N realign.hg19
#$ -o realign.hg19.out.$JOB_ID
#$ -j y
#$ -pe bscb 32
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

## In case it hasn't been previously mounted.
/programs/bin/labutils/mount_server cbsubscb08 /storage

## Set some enviroment variables.
STARTDIR=`pwd`
GENOME=/home/cgd24/nextgen/home/cgd24/genomes/hg19/bwa.rRNA-0.7.5a-r405/

## Move to the TMP directory.
SCRATCH=/SSD/cgd24_align/
mkdir $SCRATCH

## Copy data files to the working directory.
cp $GENOME/hg19* $SCRATCH
cp ~/nextgen/projects/GROseq/RawSequenceFiles/jurkat_kevin/noadapt/*.fastq.gz $SCRATCH

cd $SCRATCH

## Align reads.
for i in `ls *.fastq.gz`
do
  ## Align using BWA.
  bwa aln -t 32 hg19.rRNA $i > $i.sai
  bwa samse -n 1 -f $i.sam hg19.rRNA $i.sai $i
  samtools view -b -S $i.sam > $i.bam
  samtools sort -@ 32 $i.bam $i.sort

  ## Copy to start directory.
  cp $i.sort.bam $STARTDIR
done

rm -rf $SCRATCH
