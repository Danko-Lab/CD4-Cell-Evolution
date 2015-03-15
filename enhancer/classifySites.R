#!/usr/bin/bash

## Classifies sites as either (1) Promoter == proximal, stable.
##            :OR:            (2) Enhancer == distal (>10kb), unstable.

## Get the location of sites.
tss <- read.table("../tss_caller/HCM-U-PI.dREG-tss-clusters.tsv")

## Reposition to -110bp from the highest site.
## TODO: READ COMBINED DATA.  

## Classify sites as stable or unstable ... based on hg19.
require(stabilityHMM)

source("../lib/hmm.ss5.pa.common.R", chdir=TRUE)
seqLen= 1000

tss_hg19 <- read.table("curated_tss.bed")
dna_hg19 <- make.ss5.pA.hmm.data(collect.sequence(hg19, tss_hg19, seq.length = seqLen))
stabScore <- unstable.score(hmm= res, data= dna_hg19, n.threads = 1)

## Then classify sites as proximal/ distal based on RefGene annotations (don't want enhancer TU in the annotations).

