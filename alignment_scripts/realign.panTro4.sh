#$ -S /bin/bash
#$ -cwd
#$ -N realign.panTro4
#$ -o realign.panTro4.out.$JOB_ID
#$ -j y
#$ -pe bscb 16
#$ -M dankoc@gmail.com
#$ -m be
#$ -l h_rt=24:00:00

## In case it hasn't been previously mounted.
/programs/bin/labutils/mount_server cbsubscb08 /storage

## Set some enviroment variables.
STARTDIR=`pwd`
GENOME=/home/cgd24/nextgen/home/cgd24/genomes/panTro4/bwa.rRNA-0.7.5a-r405/

## Move to the TMP directory.
mkdir /SSD/cgd24
SCRATCH=/SSD/cgd24/$JOB_ID/
mkdir $SCRATCH

## Copy data files to the working directory.
cp $GENOME/panTro4* $SCRATCH
#cp ~/nextgen/projects/GROseq/RawSequenceFiles/CD4_nhp1_reseq/combine_pre_new_seq/C2*.fastq.gz $SCRATCH
cp ~/nextgen/projects/GROseq/RawSequenceFiles/CD4_nhp2_reseq/combine_pre_new_seq/C3*.fastq.gz $SCRATCH
#cp ~/nextgen/projects/GROseq/RawSequenceFiles/CD4_nhp3_reseq/combine_pre_new_seq/C4*.fastq.gz $SCRATCH

cd $SCRATCH

## Align reads.
for i in `ls *.fastq.gz`
do
  ## Align using BWA.
  bwa aln -t 32 panTro4.rRNA $i > $i.sai
  bwa samse -n 1 -f $i.sam panTro4.rRNA $i.sai $i
  samtools view -b -S $i.sam > $i.bam
  samtools sort -@ 32 $i.bam $i.sort

  ## Copy to start directory.
  cp $i.sort.bam $STARTDIR
done

rm -rf $SCRATCH
