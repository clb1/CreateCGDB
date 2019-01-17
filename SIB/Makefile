include ../../Makefile.common

ROOT_UCSC_URL = ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database

default:
	echo ""

download_all:
	ncftpget ${ROOT_UCSC_URL}/sibGene.txt.gz


CGDB_files.txt: SIB_hg38.gtf.gz
	echo $< > $@


browser_files.txt: SIB_hg38_Browser.gtf.gz
	echo $< > $@


SIB_hg38_unnanotated.gtf.gz: sibGene.txt.gz
	zcat $< | cut -f2- | ${UCSC_TOOLS}/genePredToGtf file stdin $(basename $@)
	gzip $(basename $@)


SIB_hg38.gtf.gz: SIB_hg38_unnanotated.gtf.gz
	zcat /raid1/projects/CGDB/models/database_merging/v2.0/GENCODE/gencode.v26.annotation.gtf.gz | egrep -v "^#" | awk '($$3=="exon")' > a_file.gtf
	zcat $< | awk '($$3=="exon")' > b_file.gtf
	intersectBed -wo -s -a a_file.gtf -b b_file.gtf > SIB_GENCODE.intersect
	scripts/annotateWithGENCODE.py SIB_GENCODE.intersect $< $@
	rm a_file.gtf b_file.gtf SIB_GENCODE.intersect


SIB_hg38_Browser.gtf.gz: SIB_hg38.gtf.gz
	echo "track type=gtf name=SIB description=SIB visibility=full" > $(basename $@)
	zcat SIB_hg38.gtf.gz >> $(basename $@)
	gzip $(basename $@)
