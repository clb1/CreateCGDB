include ../../Makefile.common

ROOT_RNACENTRAL_URL = ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/7.0/genome_coordinates

default:
	${MAKE} RNAcentral_hg38.gff3.gz


download_all:
	ncftpget ${ROOT_RNACENTRAL_URL}/Homo_sapiens.GRCh38.gff.gz


CGDB_files.txt: RNAcentral_hg38.gff.gz
	echo $< > $@


browser_files.txt: RNAcentral_hg38_Browser.gff.gz
	echo $< > $@


RNAcentral_hg38.gff.gz: Homo_sapiens.GRCh38.gff.gz
	zcat $< | egrep "^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)\b" | awk '{printf "chr%s\n",$$0;}' | gzip -c > $@


RNAcentral_hg38_Browser.gff.gz: RNAcentral_hg38.gff.gz
	echo "track name=RNAcentral description=RNAcentral visibility=full" > $(basename $@)
	zcat $< >> $(basename $@)
	gzip $(basename $@)
