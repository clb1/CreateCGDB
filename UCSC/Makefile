include ../../Makefile.common

ROOT_UCSC_URL = ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database

default:
	${MAKE} UCSC_hg38.gtf.gz


download_all:
	ncftpget ${ROOT_UCSC_URL}/knownGene.txt.gz
#	ncftpget ${ROOT_UCSC_URL}/geneName.txt.gz
#	ncftpget ${ROOT_UCSC_URL}/knownToGencodeV20.txt.gz
#	ncftpget ${ROOT_UCSC_URL}/knownToRefSeq.txt.gz
#	ncftpget ${ROOT_UCSC_URL}/knownCanonical.txt.gz Does this one have the isoform ID #'s that link to other tables?
#	ncftpget ${ROOT_UCSC_URL}/knownIsoforms.txt.gz  Does this one have the isoform ID #'s that link to other tables?
#	ncftpget ${ROOT_UCSC_URL}/knownToLocusLink.txt.gz Yes, has EntrezID
#	ncftpget ${ROOT_UCSC_URL}/refFlat.txt.gz	Yes? Is info in these three in RefSeq gtf?
#	ncftpget ${ROOT_UCSC_URL}/refGene.txt.gz	Yes?
#	ncftpget ${ROOT_UCSC_URL}/refLink.txt.gz	Yes?

CGDB_files.txt: UCSC_hg38.gtf.gz
	echo $< > $@

browser_files.txt: UCSC_hg38_Browser.gtf.gz
	echo $< > $@

UCSC_hg38.gtf.gz: knownGene.txt.gz
	zcat $< | cut -f1-10 | genePredToGtf file stdin $(basename $@)
	gzip $(basename $@)

UCSC_hg38_Browser.gtf.gz: UCSC_hg38.gtf.gz
	echo "track type=gtf name=UCSC description=UCSC visibility=full" > $(basename $@)
	zcat $< >> $(basename $@)
	gzip $(basename $@) 
