run_analysis:
	rm -f chage_expr/*
	bash makeNetGapFiles.bsh
	bash getGenesInBreaks.bsh 
	bash getCounts.bsh
	R --no-save < writeDendrogram.R 
	R --no-save < getExprChanges.limma.R
	R --no-save < analyzeTuTypes.R
	R --no-save < analyzeRearrangements.R
