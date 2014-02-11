run_analysis:
	rm chage_expr/*
	bash makeNetGapFiles.bsh
	bash getGenesInBreaks.bsh 
	bash getCounts.bsh
	R --no-save < writeDendrogram.R 
	R --no-save < getExprChanges.limmaVoom.R
	R --no-save < compareChangeFrequency.R
