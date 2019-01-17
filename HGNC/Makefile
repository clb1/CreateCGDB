#download_all:
#	wget http://www.genenames.org/cgi-bin/genefamilies/download-all/tsv

# HGNC.tsv.gz obtained by download from Custom Downloads of genenames.org

HGNC.tsv.gz: tsv
	mv $< $(basename $@)
	gzip $(basename $@)
