#!/bin/sh
#SBATCH -J chimp.unmap     # Job Name.
#SBATCH -o align.o%j      # Name of the output file (eg. myMPI.oJobID).
#SBATCH -e align.err%j    # Direct error to the error file.
#SBATCH -p normal         # Queue name.
#SBATCH -t 24:00:00       # Run time (hh:mm:ss) - 24 hours.
#SBATCH -N 1              # Requests 1 MPI node.
#SBATCH -n 3		  # 3 tasks total. ... up to 16?!
#SBATCH -A TG-MCB130131   # Siepel lab account ID to bill.
#SBATCH --mail-user=dankoc@gmail.com
#SBATCH --mail-type=all

##Set path
export PATH=$PATH:$HOME/src/genometools/genometools-1.5.1/bin/

#copy data files to the working directory
cp $WORK/genomes/panTro4/panTro4.rRNA.fa.gz /tmp/pt4.fa.gz
cd /tmp

## Index the databse.
gt suffixerator -dna -pl -tis -suf -lcp -v -parts 4 -db /tmp/pt4.fa.gz -indexname $SCRATCH/reads

## Make indexes.
gt tallymer mkindex -mersize 30 -minocc 2  -indexname $SCRATCH/30mers  -counts -pl -esa $SCRATCH/reads &
gt tallymer mkindex -mersize 50 -minocc 2  -indexname $SCRATCH/50mers  -counts -pl -esa $SCRATCH/reads &
gt tallymer mkindex -mersize 25 -minocc 2  -indexname $SCRATCH/25mers  -counts -pl -esa $SCRATCH/reads &

wait

## Write data into a useful file.
gt tallymer search -output qseqnum qpos counts sequence -strand fp -tyr $SCRATCH/30mers -q /tmp/pt4.fa.gz | gzip > $SCRATCH/30mers.gtTxt.gz &
gt tallymer search -output qseqnum qpos counts sequence -strand fp -tyr $SCRATCH/50mers -q /tmp/pt4.fa.gz | gzip > $SCRATCH/50mers.gtTxt.gz &
gt tallymer search -output qseqnum qpos counts sequence -strand fp -tyr $SCRATCH/25mers -q /tmp/pt4.fa.gz | gzip > $SCRATCH/25mers.gtTxt.gz &

wait

cp $SCRATCH/*gtTxt.gz $WORK/genomes/panTro4/unmap/

#zcat $SCRATCH/30mers.gtTxt.gz | perl ~/perl/tallymer2bed.pl | gzip > 30mer-UnMAQable.chr_coord.tsv.gz &
#zcat $SCRATCH/50mers.gtTxt.gz | perl ~/perl/tallymer2bed.pl | gzip > 50mer-UnMAQable.chr_coord.tsv.gz &
#zcat $SCRATCH/25mers.gtTxt.gz | perl ~/perl/tallymer2bed.pl | gzip > 25mer-UnMAQable.chr_coord.tsv.gz &
#
#wait

## http://www.tacc.utexas.edu/user-services/user-guides/stampede-user-guide#running

