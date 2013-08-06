run_analysis:
	bash makeNetGapFiles.bsh
	bash getCounts.bsh
	~/bin/R/R-3.0.1/bin/R --no-save < writeDendrogram.R 
	~/bin/R/R-3.0.1/bin/R --no-save < getExprChanges.R
	~/bin/R/R-3.0.1/bin/R --no-save < compareChangeFrequency.R
