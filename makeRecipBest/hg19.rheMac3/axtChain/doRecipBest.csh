#!/bin/csh -efx
# This script was automatically generated by /home/dankoc/src/kent/kent/src/hg/utils/automation/doRecipBest.pl
# It is to be executed on swiftgen.bscb.cornell.edu in /usr/projects/GROseq/NHP/makeRecipBest/hg19.rheMac3/axtChain .
# It nets in both directions to get reciprocal best chains and nets.
# This script will fail if any of its commands fail.

cd /usr/projects/GROseq/NHP/makeRecipBest/hg19.rheMac3/axtChain

# Swap hg19-best chains to be rheMac3-referenced:
chainStitchId hg19.rheMac3.over.chain.gz stdout \
| chainSwap stdin stdout \
| chainSort stdin rheMac3.hg19.tBest.chain

# Net those on rheMac3 to get rheMac3-ref'd reciprocal best net:
chainPreNet rheMac3.hg19.tBest.chain \
  ~/genomes/{rheMac3,hg19}/chrom.sizes stdout \
| chainNet -minSpace=1 -minScore=0 \
  stdin ~/genomes/{rheMac3,hg19}/chrom.sizes stdout /dev/null \
| netSyntenic stdin stdout \
| gzip -c > rheMac3.hg19.rbest.net.gz

# Extract rheMac3-ref'd reciprocal best chain:
netChainSubset rheMac3.hg19.rbest.net.gz rheMac3.hg19.tBest.chain stdout \
| chainStitchId stdin stdout \
| gzip -c > rheMac3.hg19.rbest.chain.gz

# Swap to get hg19-ref'd reciprocal best chain:
chainSwap rheMac3.hg19.rbest.chain.gz stdout \
| chainSort stdin stdout \
| gzip -c > hg19.rheMac3.rbest.chain.gz

# Net those on hg19 to get hg19-ref'd reciprocal best net:
chainPreNet hg19.rheMac3.rbest.chain.gz \
  ~/genomes/{hg19,rheMac3}/chrom.sizes stdout \
| chainNet -minSpace=1 -minScore=0 \
  stdin ~/genomes/{hg19,rheMac3}/chrom.sizes stdout /dev/null \
| netSyntenic stdin stdout \
| gzip -c > hg19.rheMac3.rbest.net.gz

# Clean up the one temp file and make md5sum:
rm rheMac3.hg19.tBest.chain
md5sum *.rbest.*.gz > md5sum.rbest.txt

# Create files for testing coverage of *.rbest.*.
netToBed -maxGap=1 rheMac3.hg19.rbest.net.gz rheMac3.hg19.rbest.net.bed
netToBed -maxGap=1 hg19.rheMac3.rbest.net.gz hg19.rheMac3.rbest.net.bed
chainToPsl rheMac3.hg19.rbest.chain.gz \
  ~/genomes/{rheMac3,hg19}/chrom.sizes \
  /gbdb/rheMac3/rheMac3.2bit /gbdb/hg19/hg19.2bit \
  rheMac3.hg19.rbest.chain.psl
chainToPsl hg19.rheMac3.rbest.chain.gz \
  ~/genomes/{hg19,rheMac3}/chrom.sizes \
  /gbdb/hg19/hg19.2bit /gbdb/rheMac3/rheMac3.2bit \
  hg19.rheMac3.rbest.chain.psl

# Verify that all coverage figures are equal:
set tChCov = `awk '{print $19;}' hg19.rheMac3.rbest.chain.psl | sed -e 's/,/\n/g' | awk 'BEGIN {N = 0;} {N += $1;} END {printf "%d\n", N;}'`
set qChCov = `awk '{print $19;}' rheMac3.hg19.rbest.chain.psl | sed -e 's/,/\n/g' | awk 'BEGIN {N = 0;} {N += $1;} END {printf "%d\n", N;}'`
set tNetCov = `awk 'BEGIN {N = 0;} {N += ($3 - $2);} END {printf "%d\n", N;}' hg19.rheMac3.rbest.net.bed`
set qNetCov = `awk 'BEGIN {N = 0;} {N += ($3 - $2);} END {printf "%d\n", N;}' rheMac3.hg19.rbest.net.bed`
if ($tChCov != $qChCov) then
  echo "Warning: hg19 rbest chain coverage $tChCov != rheMac3 $qChCov"
endif
if ($tNetCov != $qNetCov) then
  echo "Warning: hg19 rbest net coverage $tNetCov != rheMac3 $qNetCov"
endif
if ($tChCov != $tNetCov) then
  echo "Warning: hg19 rbest chain coverage $tChCov != net cov $tNetCov"
endif

mkdir experiments
mv *.bed *.psl experiments
# Make rbest net axt's download: one .axt per hg19 seq.
netSplit hg19.rheMac3.rbest.net.gz rBestNet
chainSplit rBestChain hg19.rheMac3.rbest.chain.gz
cd ..
mkdir axtRBestNet
foreach f (axtChain/rBestNet/*.net)
    netToAxt $f axtChain/rBestChain/$f:t:r.chain \
    /gbdb/hg19/hg19.2bit /gbdb/rheMac3/rheMac3.2bit stdout \
    | axtSort stdin stdout \
    | gzip -c > axtRBestNet/$f:t:r.hg19.rheMac3.net.axt.gz
  end

# Make rbest mafNet for multiz: one .maf per hg19 seq.
mkdir mafRBestNet
foreach f (axtRBestNet/*.hg19.rheMac3.net.axt.gz)
    axtToMaf -tPrefix=hg19. -qPrefix=rheMac3. $f \
        ~/genomes/{hg19,rheMac3}/chrom.sizes \
            stdout \
      | gzip -c > mafRBestNet/$f:t:r:r:r:r:r.maf.gz
end