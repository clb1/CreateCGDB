include ../../Makefile.common

SCRIPTS = /projects/CGDB2/database_merging/scripts
ROOT_ALTSCAN_URL = http://cbb.sjtu.edu.cn/~ccwei/pub/software/ALTSCAN
LIFTOVER_CHAIN_FILE = ../../hg19ToHg38.over.chain.gz 

default:
	${MAKE} ALTSCAN_hg38.gtf.gz


download_all:
	wget ${ROOT_ALTSCAN_URL}/VHC.gtf.gz
	wget ${ROOT_ALTSCAN_URL}/VMC.gtf.gz


CGDB_files.txt: ALTSCAN_hg38.gtf.gz
	echo $< > $@


browser_files.txt: ALTSCAN_hg38_Browser.gtf.gz
	echo $< > $@


ALTSCAN_hg38_Browser.gff.gz: ALTSCAN_hg38.gtf.gz
	echo "track type=gff name=ALTSCAN description=ALTSCAN visibility=full" > $(basename $@)
	zcat $< >> $(basename $@)
	gzip $(basename $@)


# Need to manually do these first:
#export PYTHONPATH=$PYTHONPATH:/usr/local/src/CrossMap/usr/lib64/python2.7/site-packages
#export PATH=$PATH:/usr/local/src/CrossMap/usr/bin
ALTSCAN_hg38_unnanotated.gtf.gz: VHC.gtf.gz VMC.gtf.gz
	zcat $^ | awk '($$3=="CDS")' | tr -d \"  > VHC_VMC.gtf
	python /usr/local/src/CrossMap/usr/bin/CrossMap.py gff ${LIFTOVER_CHAIN_FILE} VHC_VMC.gtf | grep -F fail | sed 's/.*transcript_id \(.*\).*/\1/' | \
	cut -f 1 -d';' > transcripts_failing_CrossMap.txt
	python /usr/local/src/CrossMap/usr/bin/CrossMap.py gff ${LIFTOVER_CHAIN_FILE} VHC_VMC.gtf ALTSCAN_hg38.gtf
	egrep -v "^chr[[:alnum:]]+_" ALTSCAN_hg38.gtf | scripts/removeTranscriptsFailingCrossMap.py transcripts_failing_CrossMap.txt | gzip -c > $@
	rm VHC_VMC.gtf.gz

ALTSCAN_hg38.gtf.gz: ALTSCAN_hg38_unnanotated.gtf.gz
	zcat ../RefSeq/RefSeq_hg38.gff3.gz | egrep -v "^#" | awk '{OFS="\t"}($$3=="exon"){print ($$1,$$4,$$5,$$9,"0",$$7)}' > a_file.bed
	zcat $< | awk '{OFS="\t"}($$3=="CDS"){print ($$1,$$4,$$5,$$12,"0",$$7)}' > b_file.bed
	intersectBed -wo -s -a a_file.bed -b b_file.bed > ALTSCAN_RefSeq.intersect
	scripts/annotateWithRefSeq.py ALTSCAN_RefSeq.intersect $< $@
	rm a_file.bed b_file.bed ALTSCAN_RefSeq.intersect
