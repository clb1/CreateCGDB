include ../../Makefile.common

SCRIPTS = /projects/CGDB2/database_merging/scripts
GAL_BIN = ${SCRIPTS}/GAL_0.2.2/bin
ROOT_ACEVIEW_URL = ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Aug10.human.genes
LIFTOVER_CHAIN_FILE = ../../hg19ToHg38.over.chain.gz 

default:
	${MAKE} AceView_hg38.gff.gz

download_all:
	ncftpget ${ROOT_ACEVIEW_URL}/AceView.ncbi_37.genes_gff.gff.gz
	ncftpget ${ROOT_ACEVIEW_URL}/AceView.ncbi_37.mRNA2GeneID2NM.txt.gz
	ncftpget ${ROOT_ACEVIEW_URL}/AceView.ncbi_37.good_mrnas_dna.fasta.gz
	ncftpget ${ROOT_ACEVIEW_URL}/AceView.ncbi_37.geneStructureAndSupport.txt.gz

CGDB_files.txt: AceView_hg38.gff.gz
	echo $< > $@


browser_files.txt: AceView_hg38_Browser.gff.gz
	echo $< > $@


AceView_hg38_Browser.gff.gz: AceView_hg38.gff.gz
	echo "track type=gff name=Aceview description=Aceview visibility=full" > $(basename $@)
	zcat $< >> $(basename $@)
	gzip $(basename $@)


mRNA_ends_support.tsv: AceView.ncbi_37.geneStructureAndSupport.txt.gz
	echo -e "Transcript\t5p_start_supported\t3p_stop_supported" > $@
	zcat $< |egrep "^mRNA" | awk '{OFS='\t'}{a = ($$6=="-") ? "N" : "Y"; b = ($$7=="-") ? "N" : "Y"; printf "%s\t%s\t%s\n",$$2,a,b}' >> $@


# Need to manually do these first:
#export PYTHONPATH=$PYTHONPATH:/usr/local/src/CrossMap/usr/lib64/python2.7/site-packages
#export PATH=$PATH:/usr/local/src/CrossMap/usr/bin
AceView_hg38_unannotated.gff.gz: AceView.ncbi_37.genes_gff.gff.gz mRNA_ends_support.tsv
	zcat $< | egrep -v "^mito" | awk '($$3=="exon"){printf "chr%s\n",$$0;}' | scripts/removeSuspectTranscripts.py mRNA_ends_support.tsv > AceView.ncbi_37.genes_reformatted.gff
	python /usr/local/src/CrossMap/usr/bin/CrossMap.py gff ${LIFTOVER_CHAIN_FILE} AceView.ncbi_37.genes_reformatted.gff | grep -F fail | sed 's/.*transcript_id \(.*\).*/\1/' | \
	cut -f 1 -d';' > transcripts_failing_CrossMap.txt
	python /usr/local/src/CrossMap/usr/bin/CrossMap.py gff ${LIFTOVER_CHAIN_FILE} AceView.ncbi_37.genes_reformatted.gff AceView_hg38.gff
	egrep -v "^chr[[:alnum:]]+_" AceView_hg38.gff | scripts/removeTranscriptsFailingCrossMap.py transcripts_failing_CrossMap.txt | gzip -c > $@
	rm AceView.ncbi_37.genes_reformatted.gff AceView_hg38.gff


AceView_hg38.gff.gz: AceView_hg38_unnanotated.gff.gz
	zcat /raid1/projects/CGDB/models/database_merging/v2.0/GENCODE/gencode.v26.annotation.gtf.gz | egrep -v "^#" | awk '($$3=="exon")' > a_file.gtf
	zcat $< | awk '($$3=="exon")' > b_file.gtf
	intersectBed -wo -s -a a_file.gtf -b b_file.gtf > AceView_GENCODE.intersect
	scripts/annotateWithGENCODE.py AceView_GENCODE.intersect $< $@
	rm a_file.gtf b_file.gtf AceView_GENCODE.intersect


#AceView_hg38.gff3.gz: AceView_hg38.gff.gz AceView.ncbi_37.mRNA2GeneID2NM.txt.gz
#	${GAL_BIN}/gtf2gff3 <(zcat $<) | ${SCRIPTS}/addAliasesToAceviewGFF3.py AceView.ncbi_37.mRNA2GeneID2NM.txt.gz $@
