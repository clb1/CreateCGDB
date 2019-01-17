#!/usr/bin/env python

from collections import defaultdict
import gzip
import re
import sys


def parseBiomartTable(transcript_aliases, biomart_tsv, labels):
    with gzip.open(biomart_tsv, "rb") as ip:
        header = ip.readline()
        for line in ip:
            tup = line.strip("\n").split("\t")
            for i in xrange(len(labels)):
                if (tup[i+2] != ""):
                    transcript_aliases[tup[1]].add( "%s%s" % (labels[i], tup[i+2]) )


def mapTranscriptToGene(gencode_transcript_to_gene):
    transcript_to_gene = {}
    with open(gencode_transcript_to_gene, 'r') as ip:
        for line in ip:
            gencode_transcript, gencode_gene = line.strip().split()
            transcript_to_gene[gencode_transcript] = gencode_gene
    return transcript_to_gene

    
def getGeneEntrezAlias(transcript_to_gene, transcript_to_entrez, gene_aliases):
    with gzip.open(transcript_to_entrez, 'rb') as ip:
        for line in ip:
            gencode_transcript, entrezID = line.strip().split("\t")
            if (transcript_to_gene.has_key(gencode_transcript)):
                gene_aliases[ transcript_to_gene[gencode_transcript] ].add(entrezID)


def getGeneHGNCAlias(transcript_to_gene, transcript_to_hgnc, gene_aliases):
    with gzip.open(transcript_to_hgnc, 'rb') as ip:
        for line in ip:
            gencode_transcript, gene_symbol = line.strip().split("\t")
            if (transcript_to_gene.has_key(gencode_transcript)):
                gene_aliases[ transcript_to_gene[gencode_transcript] ].add(gene_symbol)


def getTranscriptRefseqAlias(transcript_to_refseq, transcript_aliases):
    with gzip.open(transcript_to_refseq, 'rb') as ip:
        for line in ip:
            tup = line.strip().split("\t")
            transcript_aliases[ tup[0] ].add(tup[1])


def modifyAndWrite(transcript_aliases, gene_aliases, input_gff3, output_gff3):
    ip = gzip.open(input_gff3, 'rb')
    op = gzip.open(output_gff3, 'wb')

    for line in ip:
        fields = line.strip().split("\t")
        if (line[0] == '#'):
            op.write(line)
        #elif (fields[2] == "gene"):
            #gencode_gene = fields[-1][3:].split(';')[0].split('.')[0]
            #if ("gene_name=" in fields[-1]):
            #    gene_name = fields[-1].split("gene_name=")[1].split(';')[0]
            #    assert (len(gene_name) > 0)
            #    gene_aliases[gencode_gene].add(gene_name)
            #if (gene_aliases.has_key(gencode_gene)):
            #    #fields[-1] += ";Alias=%s" % ",".join(list(gene_aliases[gencode_gene]))
            #    fields[-1] = ";Alias=%s" % gene_aliases[gencode_gene]
            #    op.write("%s\n" % "\t".join(fields))
        elif (fields[2] == "transcript"):
            gencode_transcript = fields[-1][3:].split(';')[0].split('.')[0]
            if (transcript_aliases.has_key(gencode_transcript)):
                fields[-1] += ";Alias=%s" % ",".join(list(transcript_aliases[gencode_transcript]))
                op.write("%s\n" % "\t".join(fields))
        else:
            op.write(line)
                
    ip.close()
    op.close()
    

if (__name__ == "__main__"):
    biomart_Entrez_HGNC_tsv, biomart_RNACentral_UCSC_tsv, biomart_RefSeq_tsv, input_gff3, output_gff3 = sys.argv[1:]
    
    gene_aliases = {}
    transcript_aliases = defaultdict(set)
    for biomart_tsv, labels in [(biomart_Entrez_HGNC_tsv, ["Entrez:",""]),
                                (biomart_RNACentral_UCSC_tsv, ["RNACentral:", "UCSC:"]),
                                (biomart_RefSeq_tsv, ["RefSeq:", "RefSeq:"])]:
            parseBiomartTable(transcript_aliases, biomart_tsv, labels)

    transcript_to_gene = mapTranscriptToGene(gencode_transcript_to_gene)
    getGeneEntrezAlias(transcript_to_gene, transcript_to_entrez, gene_aliases)
    getGeneHGNCAlias(transcript_to_gene, transcript_to_hgnc, gene_aliases)
    getTranscriptRefseqAlias(transcript_to_refseq, transcript_aliases)

    modifyAndWrite(transcript_aliases, gene_aliases, input_gff3, output_gff3)
    
    sys.exit(0)
    
