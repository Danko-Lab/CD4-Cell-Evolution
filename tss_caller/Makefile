run_1:
	#bash nhptss_liftOver.bsh
	#bash testSpec.bsh
	bash bedGraphToBed.dreghd.bsh
	bash prepStabilityClassification.bsh
	R --no-save < classifySites.R

run_2:
	bash getHumanGains.bsh
	bash runInsight.bsh
	bash getMotifEnrichment.bsh

