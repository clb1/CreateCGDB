include ../../Makefile.common

ROOT_FANTOM5_URL = http://fantom.gsc.riken.jp/5/datafiles

default:
	${MAKE} hg38.cage_peak_coord_robust.bed.gz
	${MAKE} hg38.cage_peak_coord_permissive.bed.gz
#	${MAKE} hg38.cage_peak_phase1and2combined_coord.bed.gz


download_all:
	wget ${ROOT_FANTOM5_URL}/phase1.3/extra/CAGE_peaks/hg19.cage_peak_coord_permissive.bed.gz
	wget ${ROOT_FANTOM5_URL}/phase1.3/extra/CAGE_peaks/hg19.cage_peak_coord_robust.bed.gz
#	wget ${ROOT_FANTOM5_URL}/phase2.0/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz


# hg38.cage_peak_phase1and2combined_coord.bed.gz
# hg38.cage_peak_coord_permissive.bed.gz
CGDB_files.txt: hg38.cage_peak_coord_robust.bed.gz
	-rm $@
	for f in $^; do \
		echo $$f >> $@; \
	done


# 		hg38.cage_peak_phase1and2combined_coord_Browser.bed.gz
browser_files.txt: hg38.cage_peak_coord_robust_Browser.bed.gz hg38.cage_peak_coord_permissive_Browser.bed.gz
	-rm $@
	for f in $^; do \
		echo $$f >> $@; \
	done


hg38.cage_peak_%.bed.gz: hg19.cage_peak_%.bed.gz
	${CROSSMAP} bed ${LIFTOVER_CHAIN_FILE} <(zcat $<) hg38_$*.bed
	egrep -v "alt|random|chrUn" hg38_$*.bed | cut -f 1-6 | sort -k1,1 -k2,2n | \
	bedtools merge -s -c 4,5 -o collapse,sum -i stdin | awk '{OFS="\t"}{print $$1,$$2,$$3,$$5,$$6,$$4}' | gzip -c > $@
	rm hg38_$*.bed


hg38.cage_peak_%_Browser.bed.gz: hg38.cage_peak_%.bed.gz
	echo "track type=bedDetail name=FANTOM_$* description=FANTOM_$* visibility=full" > $(basename $@)
	zcat $< >> $(basename $@)
	gzip $(basename $@)
