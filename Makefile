SHELL = /bin/bash

# VERSION is to be defined by an including Makefile
VERSION = $(shell c=`pwd`; basename $$c)

ROOT_DIR = /raid1/projects/CGDB/models/database_merging
SCRIPTS = ${ROOT_DIR}/scripts
COMMON_BIN = ${ROOT_DIR}/bin

# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# Use for CrossMap (http://crossmap.sourceforge.net/) and liftOver
LIFTOVER_CHAIN_FILE = ${ROOT_DIR}/hg19ToHg38.over.chain.gz

CROSSMAP = export PYTHONPATH=/usr/local/src/CrossMap/usr/lib64/python2.7/site-packages; python /usr/local/src/CrossMap/usr/bin/CrossMap.py

JBROWSE_ROOT_DIR = /srv/www/htdocs/jbrowse/JBrowse-1.11.6
JBROWSE_INSTALL_DIR = ${JBROWSE_ROOT_DIR}/sample_data/json/human

ROOT_URL = http://ism.ucsd.edu/CGDB_Sources

INDEXED_GENOME_FASTA = /raid1/references_and_indexes/hg38/hg38.fa
TSS_SOURCES = FANTOM5 # DBTSS
ISOFORM_MODEL_SOURCES = SIB AceView HInv ALTSCAN GENCODE RefSeq RNAcentral UCSC
POLYA_SOURCES = APADB # polyAMotif xPAD  # polyADB2 APASdb PolyA-seq

ALL_SOURCES = $(TSS_SOURCES) $(ISOFORM_MODEL_SOURCES) $(POLYA_SOURCES)

TEMP_DIR = /raid1/projects/scratch

.PHONY: CGDB convert_all_gff_to_bed download_all browser_files.txt


CGDB:
    ${MAKE} CGDB_${VERSION}.bed


# For STAR. Here because main CGDB GFF hasn't been created yet. Maybe move to main CGDB directory's Makefile once the CGDB GFF3 is created.
all_SJs.tsv.bz2:
    i=''; \
    for s in $(ISOFORM_MODEL_SOURCES); do \
        ${MAKE} -C $$s CGDB_files.txt; \
        for f in `cat $$s/CGDB_files.txt`; do \
            i="$$i $$s/$$f"; \
        done; \
    done; \
    i=$(strip $$i); \
    i=`echo $$i | sed 's/ /,/g'`; \
    ${SCRIPTS}/compileAllSJs.py $$i $@


# For creating microexon splice junctions for use with STAR. Same comment as above.
all_microexon_SJs.tsv: microexons/microexons_hg38.bed
    i=''; \
    for s in $(ISOFORM_MODEL_SOURCES); do \
        ${MAKE} -C $$s CGDB_files.txt; \
        for f in `cat $$s/CGDB_files.txt`; do \
            i="$$i $$s/$$f"; \
        done; \
    done; \
    i=$(strip $$i); \
    i=`echo $$i | sed 's/ /,/g'`; \
    ${SCRIPTS}/createMicroexonsSJsTSV.py $< $$i $@


#CGDB_${VERSION}.bed CGDB${VERSION}.gtf CGDB${VERSION}.gff3 CGDB${VERSION}.pkl CGDB${VERSION}_absorption_info.txt
CGDB: hg38_chromosomes_sizes.tsv HGNC/HGNC.tsv.gz
    @t=''; \
    i=''; \
    p=''; \
    for s in $(TSS_SOURCES); do \
        ${MAKE} -C $$s CGDB_files.txt; \
        for f in `cat $$s/CGDB_files.txt`; do \
            t="$$t $$s/$$f"; \
        done; \
    done; \
    t=$(strip $$t); \
    t=`echo $$t | sed 's/ /,/g'`; \
    for s in $(ISOFORM_MODEL_SOURCES); do \
        ${MAKE} -C $$s CGDB_files.txt; \
        for f in `cat $$s/CGDB_files.txt`; do \
            i="$$i $$s/$$f"; \
        done; \
    done; \
    i=$(strip $$i); \
    i=`echo $$i | sed 's/ /,/g'`; \
    for s in $(POLYA_SOURCES); do \
        ${MAKE} -C $$s CGDB_files.txt; \
        for f in `cat $$s/CGDB_files.txt`; do \
            p="$$p $$s/$$f"; \
        done; \
    done; \
    p=$(strip $$p); \
    p=`echo $$p | sed 's/ /,/g'`; \
    echo "scripts/createCGDB.py ${TEMP_DIR} ${INDEXED_GENOME_FASTA} $$t $$i $$p $^ CGDB${VERSION}    " ; \
    ${SCRIPTS}/createCGDB.py ${TEMP_DIR} ${INDEXED_GENOME_FASTA} $$t $$i $$p $^ CGDB${VERSION}


hg38_chromosomes_sizes.tsv:
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "select chrom, size from hg38.chromInfo"  > $@

#
###### Isommune JBrowse-related ######
#
# For the GFF3 isoform models
convert_all_gff_to_bed:
    i=''; \
    for s in $(ISOFORM_MODEL_SOURCES); do \
        ${MAKE} -C $$s CGDB_files.txt; \
        for f in `cat $$s/CGDB_files.txt`; do \
            i="$$i $$s/$$f"; \
        done; \
    done; \
    i=$(strip $$i); \
    i=`echo $$i | sed 's/ /,/g'`; \
    ${SCRIPTS}/convertToBed12.py HGNC/HGNC.tsv.gz $$i


#
###### UCSC Browser-related ######
#

CGDB2_Browser.bed: CGDB2.bed
    echo "track type=bed name=CGDB2 description=CGDB2 itemRgb=On visibility=full" > $@
    cat $< >> $@

CGDB2_Browser.gtf: CGDB2.gtf
    echo "track type=gtf name=CGDB2 description=CGDB2 useScore=1 visibility=full" > $@
    cat $< >> $@

CGDB2_sources_for_Brower.txt:
    for s in $(ALL_SOURCES); do \
        ${MAKE} -C $$s browser_files.txt; \
        for f in `cat $$s/browser_files.txt`; do \
            echo "${ROOT_URL}/$$s/$$f" >> $@ ; \
        done; \
    done

clean-browser-related-files:
    rm */browser_files.txt */*Browser*
