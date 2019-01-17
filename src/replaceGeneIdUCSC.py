#!/usr/bin/env python

import gzip
import re
import sys


def mapUCSCToGencodeGene(biomart_HGNC_Entrez_UCSC_tsv, ucsc_to_gencode_tsv):
    ucsc_to_gencode_gene = {}
    gencode_transcript_to_gene = {}
    gene_aliases = {}
    
    with gzip.open(biomart_HGNC_Entrez_UCSC_tsv, "rb") as ip:
        header = ip.readline()
        for line in ip:
            aliases = []

            gencode_gene, gencode_transcript, HGNC, entrezID, ucscID = line.strip("\n").split("\t")
            ucscID = ucscID.split('.')[0]
            gencode_transcript_to_gene[gencode_transcript] = gencode_gene

            if (HGNC != ''):
                aliases.append(HGNC)
            if (entrezID != ''):
                aliases.append( "Entrez:%s" % entrezID )

            assert (not (ucscID != '' and gencode_gene == ''))
            if (len(aliases) > 0):
                alias_string = ",".join(aliases)
                if (ucscID != ''):
                    gene_aliases[ucscID] = alias_string
                #gene_aliases[gencode_gene] = alias_string

            if (ucscID != ''):
                ucsc_to_gencode_gene[ucscID] = gencode_gene

    with gzip.open(ucsc_to_gencode_tsv, "rb") as ip:
        for line in ip:
            ucscID, gencode_transcript = map(lambda x: x.split('.')[0], line.split("\t"))
            if (ucsc_to_gencode_gene.has_key(ucscID)):
                try:
                    assert (ucsc_to_gencode_gene[ucscID] == gencode_transcript_to_gene[gencode_transcript])
                except AssertionError:
                    import pdb
                    pdb.set_trace()
                        
            elif (gencode_transcript_to_gene.has_key(gencode_transcript)):
                ucsc_to_gencode_gene[ucscID] = gencode_transcript_to_gene[gencode_transcript]

    return ucsc_to_gencode_gene, gene_aliases


def relabelAndWrite(ucsc_to_gencode_gene, gene_aliases, input_gtf, output_gtf):
    reUCSCLine = re.compile("(.+gene_id .)(uc\d\d\d\S\S\S)\.\d+(.. transcript_id .+)")

    op = gzip.open(output_gtf, 'wb')
    with gzip.open(input_gtf, 'rb') as ip:
        for line in ip:
            mo = reUCSCLine.match(line)
            if (mo):
                pre, ucscID, post = mo.groups()

                alias_string = ''
                if (gene_aliases.has_key(ucscID)):
                    alias_string = " Alias \"%s\";" % gene_aliases[ucscID]

                if (ucsc_to_gencode_gene.has_key(ucscID)):
                    gencode_gene = ucsc_to_gencode_gene[ucscID]
                    op.write("%s%s%s%s\n" % (pre, gencode_gene, post.strip(), alias_string))
                else:
                    op.write(line)
            else:
                op.write(line)
    op.close()

    
if (__name__ == "__main__"):
    ucsc_to_gencode_tsv, biomart_HGNC_Entrez_UCSC_tsv, input_gtf, output_gtf = sys.argv[1:] # ucsc_to_entrez_tsv, 

    ucsc_to_gencode_gene, gene_aliases = mapUCSCToGencodeGene(biomart_HGNC_Entrez_UCSC_tsv, ucsc_to_gencode_tsv)
    relabelAndWrite(ucsc_to_gencode_gene, gene_aliases, input_gtf, output_gtf)

    sys.exit(0)
