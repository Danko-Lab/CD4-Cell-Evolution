run_all:
	make -C Alignments_1stPrep/
	make -C ../Alignments_2ndPrep/
	make -C ../Alignments_3rdPrep/
	make -C ../AllData/
	make run_1 -C ../tss_caller/
	make -C ../tu_caller/
	make -C ../annotations/
	make run_2 -C ../tss_caller/
