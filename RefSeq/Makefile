include ../../Makefile.common

# RefSeq data is at ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens. 
# Current GFFs file is in directory GFF. Older ones are under ARCHIVE directory.
# See also http://www.ncbi.nlm.nih.gov/refseq/
#

.PHONY:

default:
	${MAKE} RefSeq_hg38.gff3.gz


interim_GRCh38.p10_top_level_2017-05-03.gff3.gz:
	ncftpget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF_interim/interim_GRCh38.p10_top_level_2017-05-03.gff3.gz
	ncftpget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/chr_accessions_GRCh38.p7
	#ncftpget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/unlocalized_accessions_GRCh38.p2
	#ncftpget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/alts_accessions_GRCh38.p2
	#ncftpget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/unplaced_accessions_GRCh38.p2


#unlocalized_accessions_GRCh38.p2 alts_accessions_GRCh38.p2 unplaced_accessions_GRCh38.p2
RefSeq_hg38.gff3.gz: interim_GRCh38.p10_top_level_2017-05-03.gff3.gz chr_accessions_GRCh38.p7
	${SCRIPTS}/convertRefSeqChromosomeNames.py $@ $^


CGDB_files.txt: RefSeq_hg38.gff3.gz
	echo $< > $@

browser_files.txt: RefSeq_hg38_Browser.gff3.gz
	echo $< > $@


RefSeq_hg38_Browser.gff3.gz: RefSeq_hg38.gff3.gz
	echo "track name=RefSeq description=RefSeq visibility=full" > $(basename $@)
	zcat $< | egrep -v "^#" >> $(basename $@)
	gzip $(basename $@)
