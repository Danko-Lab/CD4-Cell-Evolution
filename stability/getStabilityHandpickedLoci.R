#liftOver curated_tss.bed ../makeRecipBest/hg19.panTro4/axtChain/hg19.panTro4.rbest.chain.gz curated_tss.panTro4.bed unmap
#liftOver curated_tss.bed ../makeRecipBest/hg19.rheMac3/axtChain/hg19.rheMac3.rbest.chain.gz curated_tss.rheMac3.bed unmap

source("hmm.ss5.pa.common.R", chdir=TRUE)
seqLen= 1500

tss_hg19 <- read.table("curated_tss.bed")
dna_hg19 <- make.ss5.pA.hmm.data(collect.sequence(hg19, tss_hg19, seq.length = seqLen))
unstable.score(hmm= res, data= dna_hg19, n.threads = 1)

tss_panTro4 <- read.table("curated_tss.panTro4.bed")
dna_panTro4 <- make.ss5.pA.hmm.data(collect.sequence(panTro4, tss_panTro4, seq.length = seqLen))
unstable.score(hmm= res, data= dna_panTro4, n.threads = 1)

tss_rheMac3 <- read.table("curated_tss.rheMac3.bed")
dna_rheMac3 <- make.ss5.pA.hmm.data(collect.sequence(rheMac3, tss_rheMac3, seq.length = seqLen))
unstable.score(hmm= res, data= dna_rheMac3, n.threads = 1)


## Sanity check...
#[cgd24@cbsubscb08 stability]$ cat ../annotations/gencode.v18.transcript.tsv | grep "protein_coding" | grep "+" | head -n 100 > tmp.plus
summary(unstable.score(hmm= res, data= make.ss5.pA.hmm.data(collect.sequence(hg19, read.table("tmp.plus"))), n.threads = 1))
#[cgd24@cbsubscb08 stability]$ cat ../annotations/gencode.v18.transcript.tsv | grep "protein_coding" | grep "-" | head -n 100 > tmp.minus
summary(unstable.score(hmm= res, data= make.ss5.pA.hmm.data(collect.sequence(hg19, read.table("tmp.minus"))), n.threads = 1))

