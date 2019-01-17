include ../../Makefile.common

default:
	${MAKE} HInv_hg38.gff.gz

download_all:
	wget http://www.h-invitational.jp/hinv/download/map/gff3/h-inv_pub.gff3.gz

CGDB_files.txt: HInv_hg38.gff3.gz
	echo $< > $@

browser_files.txt: HInv_hg38_Browser.gff3.gz
	echo $< > $@

HInv_hg38_unnanotated.gff3.gz: h-inv_pub.gff3.gz
	zcat h-inv_pub.gff3.gz | awk '($$3=="mRNA")' | sed 's/Parent=/\t/' | cut -f 1,10 | uniq > referenced_parent_genes.txt
	${CROSSMAP} gff ${LIFTOVER_CHAIN_FILE} h-inv_pub.gff3.gz | grep fail | sed 's/;/\t/g' | cut -f 3,9 | sed 's/ID=//' > failed_CrossMap_elements.txt
	${CROSSMAP} gff ${LIFTOVER_CHAIN_FILE} <(zcat $<) HInv_hg38.gff3
	egrep -v "^chr[[:alnum:]]+_" HInv_hg38.gff3 | scripts/removeElementsFailingCrossMap.py failed_CrossMap_elements.txt referenced_parent_genes.txt | gzip -c > $@
	rm HInv_hg38.gff3

HInv_hg38.gff3.gz: HInv_hg38_unnanotated.gff3.gz
	#zcat /raid1/projects/CGDB/models/database_merging/v2.0/GENCODE/gencode.v26.annotation.gtf.gz | egrep -v "^#" | awk '($$3=="exon")' > a_file.gtf
	#zcat $< | egrep -v "^#" | awk '($$3=="exon")' > b_file.gtf
	#intersectBed -wo -s -a a_file.gtf -b b_file.gtf > HInv_GENCODE.intersect
	scripts/annotateWithGENCODE.py HInv_GENCODE.intersect $< $@
	#rm a_file.gtf b_file.gtf HInv_GENCODE.intersect

HInv_hg38_Browser.gff3.gz: HInv_hg38.gff3.gz
	-rm $@
	echo "track name=H-Inv description=H-Inv visibility=full" > $(basename $@)
	zcat $< | awk '($$3!="gene" && $$3!="mRNA")' >> $(basename $@)
	gzip $(basename $@)
