#!/usr/bin/env python

from collections import defaultdict
from operator import itemgetter
import sys


def groupCGDBIsoformLabelsByLocus(CGDB_bed):
    isoform_ID_parts = defaultdict(list)

    with open(CGDB_bed, 'r') as ip:
        for line in ip:
            fields = line.split("\t")
            if (fields[3].startswith("CG_")):
                locus, tss, splice, polyA = fields[3].split('.')
                isoform_ID_parts[locus].append( (tss, splice, polyA, fields[0], int(fields[1]), int(fields[2])) )

    return isoform_ID_parts


def writeStats(isoform_ID_parts, output_tsv):
    all_loci_data = []
    for locusID, parts in isoform_ID_parts.items():
        num_tss = len(set(map(itemgetter(0), parts)))
        num_splice = len(set(map(itemgetter(1), parts)))
        num_polyA = len(set(map(itemgetter(2), parts)))
        locus_start = min(map(itemgetter(4), parts))
        locus_stop = min(map(itemgetter(5), parts))

        chrom = set(map(itemgetter(3), parts))
        assert (len(chrom) == 1)
        chrom = list(chrom)[0]

        tss_splice_ratio = float(num_tss)/float(num_splice)
        polyA_splice_ratio = float(num_polyA)/float(num_splice)
        all_loci_data.append( (locusID, num_tss, num_splice, num_polyA, tss_splice_ratio, polyA_splice_ratio, chrom, locus_start, locus_stop) )

    all_loci_data.sort(key=itemgetter(4,5), reverse=True)
    all_loci_data.sort(key=itemgetter(2))


    with open(output_tsv, 'w') as op:
        for tup in all_loci_data:
            op.write("%s\t%s\t%s\t%s\t%4.3f\t%4.3f\t%s:%d-%d\n" % tup)


if (__name__ == "__main__"):
    CGDB_bed, output_tsv = sys.argv[1:]

    isoform_ID_parts = groupCGDBIsoformLabelsByLocus(CGDB_bed)
    writeStats(isoform_ID_parts, output_tsv)
    
    sys.exit(0)
