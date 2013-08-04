run_analysis:
	bash makeNetGapFiles.bsh
	bash getCounts.bsh
	R --no-save < writeDendrogram.R 
