#!/usr/bin/env python

import gzip
import sys


def readConversionFiles(chromosome_accessions):
    accession_to_chrom = {}
    ip = open(chromosome_accessions, 'r')
    for line in ip:
        if (line[0] != '#'):
            fields = line.strip().split("\t")
            if (fields[0] == "MT"):
                chrom = "chrM"
            else:
                chrom = "chr%s" % fields[0]
            accession_to_chrom[fields[1]] = chrom
    ip.close()
    return accession_to_chrom


def convertGFF(input_gff, accession_to_chrom, output_gff):
    ip = gzip.open(input_gff, 'rb')
    op = gzip.open(output_gff, 'wb')
    for line in ip:
        if (line.startswith("##sequence-region")):
            fields = line.split()
            if (accession_to_chrom.has_key(fields[1])):
                fields[1] = accession_to_chrom[fields[1]]
                op.write("%s\n" % " ".join(fields))
        elif (line[0] == '#'):
            op.write(line)
        else:
            accession, rest_of_line = line.split("\t", 1)
            if (accession_to_chrom.has_key(accession)):
                chrom_name = accession_to_chrom[accession]
                op.write("%s\t%s" % (chrom_name, rest_of_line))
    ip.close()
    op.close()
    

if (__name__ == "__main__"):
    output_gff, input_gff, chromosome_accessions = sys.argv[1:]

    accession_to_chrom = readConversionFiles(chromosome_accessions)
    convertGFF(input_gff, accession_to_chrom, output_gff)
    
    sys.exit(0)
    
