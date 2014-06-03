run_all:
	make -C Alignments_1stPrep/
	make -C ../Alignments_2ndPrep/
	make -C ../Alignments_3rdPrep/
	make -C ../AllData/
	#make -C ../tss_caller/ ## BOTH OF THESE REQUIRE EACH OTHER.
	#make -C ../annotations/
