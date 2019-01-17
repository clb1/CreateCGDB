#!/usr/bin/env python

import sys
import gzip


def getConversionDict(source):
    D = None
    if (source == "RefSeq"):
        D = {"CDS": "CDS", "gene": "gene", "mRNA": "mRNA",
             "ncRNA": "mRNA", "rRNA": "mRNA", "tRNA": "mRNA", "transcript": "mRNA", "primary_transcript": "mRNA",
             "exon": "exon", "match": "exon", "cDNA_match": "exon",
             "long_terminal_repeat": "mRNA", "D_loop": "gene",
             "C_gene_segment": "mRNA", "D_gene_segment": "mRNA", "J_gene_segment": "mRNA", "V_gene_segment": "mRNA"}
    else:
        print >> sys.stderr, "ERROR: unrecognized isoform model source -> %s. Exiting." % source
        sys.exit(1)

    return D

def standardize(gzip_gff3, D):
    ip = gzip.open(gzip_gff3, 'rb')
    op = sys.stdout
    for line in ip:
        if (line[0] == '#'):
            op.write(line)
        else:
            fields = line.split("\t")
            if (D.has_key(fields[2])):
                fields[2] = D[ fields[2] ]
                op.write("\t".join(fields))
            else:
                op.write(line)
    op.close()
    ip.close()
    

if (__name__ == "__main__"):
    source, gzip_gff3 = sys.argv[1:]

    D = getConversionDict(source)
    standardize(gzip_gff3, D)

    sys.exit(0)
