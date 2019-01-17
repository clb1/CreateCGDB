include ../../Makefile.common

# The main GENCODE site is http://www.gencodegenes.org.

SCRIPTS = /raid1/projects/CGDB/database_merging/scripts
RELEASE_NUMBER = 26
ROOT_GENCODE_URL = ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${RELEASE_NUMBER}

download_all:
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.annotation.gff3.gz
	ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.annotation.gtf.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.polyAs.gff3.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.2wayconspseudos.gtf.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.long_noncoding_RNAs.gtf.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.tRNAs.gtf.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.metadata.EntrezGene.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.metadata.HGNC.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.metadata.PolyA_feature.gz
	#ncftpget ${ROOT_GENCODE_URL}/gencode.v${RELEASE_NUMBER}.metadata.RefSeq.gz


GENCODE_hg38.gff3.gz: gencode.v${RELEASE_NUMBER}.annotation.gff3.gz
	ln -s $< $@


CGDB_files.txt: GENCODE_hg38.gff3.gz
	echo $< > $@

browser_files.txt: GENCODE_hg38_Browser.gff3.gz
	echo $< > $@

transcript_to_gene.tsv.gz: gencode.v${RELEASE_NUMBER}.annotation.gff3.gz
	zcat $< | awk '($$3=="transcript"){print $$9;}' | tr ';=' ' '| cut -f 2,4 -d$$' ' | gzip -c > $@

#GENCODE_hg38.gff3.gz: gencode.v${RELEASE_NUMBER}.annotation.gff3.gz gencode.v${RELEASE_NUMBER}.metadata.EntrezGene.gz \
#			gencode.v${RELEASE_NUMBER}.metadata.HGNC.gz gencode.v${RELEASE_NUMBER}.metadata.RefSeq.gz
#	${SCRIPTS}/addAliasesToGencodeGFF3.py <(zcat $< | awk '($$3=="transcript"){print $$9;}' | tr ';=' ' '| cut -f 2,4 -d$$' ') \
#	gencode.v${RELEASE_NUMBER}.metadata.EntrezGene.gz gencode.v${RELEASE_NUMBER}.metadata.HGNC.gz gencode.v${RELEASE_NUMBER}.metadata.RefSeq.gz $< $@


#GENCODE_hg38.gff3.gz: biomart_Entrez_HGNC.tsv.gz biomart_RNACentral_UCSC.tsv.gz biomart_RefSeq.tsv.gz gencode.v${RELEASE_NUMBER}.annotation.gff3.gz
#	${SCRIPTS}/addAliasesToGencodeGFF3.py $^ $@


GENCODE_hg38_Browser.gtf.gz: GENCODE_hg38.gff3.gz
	echo "track type=gff name=GENCODEv${RELEASE_NUMBER} description=GENCODEv${RELEASE_NUMBER} visibility=full" > $(basename $@)
	zcat $< | egrep -v "^#" >> $(basename $@)
	gzip $(basename $@)
