run_analysis:
#	bash makeNetGapFiles.bsh
#	bash getGenesInBreaks.bsh 
#	bash getCounts.bsh
	R --no-save < writeDendrogram.R 
	R --no-save < getExprChanges.DESeq2.R
	R --no-save < compareChangeFrequency.R
