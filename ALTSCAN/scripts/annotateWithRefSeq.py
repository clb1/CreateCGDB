#!/usr/bin/env python

from collections import defaultdict, Counter
from operator import itemgetter
import gzip
import re
import sys

import pdb

def compileALTSCANAnnotation(intersect_file):
    annot_per_transcript = defaultdict(list)
    accum_annot = defaultdict(set)

    ip = open(intersect_file, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        chrom = fields[0]
        exon_start = fields[1]
        exon_stop = fields[2]
        altscan_transcript = fields[9].replace(';','')

        attributes_dict = dict(map(lambda x: tuple(x.split('=')), fields[3].split(';')))
        for Db, Db_acc in map(lambda y: tuple(y.split(':', 1)), attributes_dict['Dbxref'].split(',')):
            if (Db == "HGNC" and not Db_acc.endswith("HGNC")): # Some entries are corrupted with "HGNC:HGNC:HGNC" in the RefSeq file
                accum_annot[altscan_transcript].add( (chrom, exon_start, exon_stop, Db, Db_acc) )

        if (attributes_dict.has_key("gene")):
            RefSeq_acc = attributes_dict["gene"]
            accum_annot[altscan_transcript].add( (chrom, exon_start, exon_stop, "RefSeq", RefSeq_acc) )

    ip.close()

    for altscan_transcript, exon_annots in accum_annot.items():
        C = Counter(map(itemgetter(3,4), exon_annots))
        most_common = C.most_common()
        if (len(most_common) == 0):
            pdb.set_trace()
        while (len(most_common) > 0):
            if (most_common[0][0][0]=="HGNC"):
                annot_per_transcript[altscan_transcript].append( most_common[0][0] )
                most_common = filter(lambda x: x[0][0]!="HGNC", most_common)
            elif (most_common[0][0][0]=="RefSeq"):
                annot_per_transcript[altscan_transcript].append( most_common[0][0] )
                most_common = filter(lambda x: x[0][0]!="RefSeq", most_common)

    return annot_per_transcript


def addAnnotaton(annot_per_transcript, unannot_gtf_gz, output_gtf_gz):
    reGTF = re.compile("^chr.+\s+altscan\s+CDS\s+\d+\s+\d+\s+.+\s+.\s+\d\s+gene_id \S+; transcript_id (\S+);")

    ip = gzip.open(unannot_gtf_gz, 'rb')
    op = gzip.open(output_gtf_gz, 'wb')

    for line in ip:
        line = line.strip()
        mo = reGTF.match(line)
        altscan_transcript = mo.group(1)

        #assert (len(annot_per_transcript[altscan_locus]) > 0), "No annotation for locus %s" % altscan_locus

        for Db, Db_acc in annot_per_transcript[altscan_transcript]:
            line += " %s %s;" % (Db, Db_acc)

        op.write("%s\n" % line)

    ip.close()
    op.close()


if (__name__ == "__main__"):
    intersect_file, unannot_gtf_gz, output_gtf_gz = sys.argv[1:]

    annot_per_transcript = compileALTSCANAnnotation(intersect_file)
    addAnnotaton(annot_per_transcript, unannot_gtf_gz, output_gtf_gz)
    sys.exit(0)
