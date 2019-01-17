#!/usr/bin/env python

from collections import defaultdict, Counter
from operator import itemgetter
import gzip
import re
import sys

import pdb

def compileAnnotation(intersect_file):
    annot_per_transcript = defaultdict(dict)
    accum_annot = defaultdict(set)

    ip = open(intersect_file, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        chrom = fields[0]
        gencode_exon_start = fields[3]
        gencode_exon_stop = fields[4]
        gencode_annot = fields[8]
        
        hinv_exon_start = fields[12]
        hinv_exon_stop = fields[13]
        hinv_annot = fields[17]
        
        if (gencode_exon_start==hinv_exon_start or gencode_exon_stop==hinv_exon_stop):
            gencode_annot_elems = gencode_annot.split()
            gencode_annot_dict = dict([(gencode_annot_elems[i],gencode_annot_elems[i+1]) for i in xrange(0,len(gencode_annot_elems),2)])
            gencode_gene = gencode_annot_dict['gene_id'][1:-2]

            hinv_transcript = hinv_annot.split()[0][7:]

            accum_annot[hinv_transcript].add( ("GENCODE", gencode_gene) )
            if (gencode_annot_dict.has_key("gene_name")):
                accum_annot[hinv_transcript].add( ("RefSeq", gencode_annot_dict["gene_name"][1:-2]) )

    ip.close()

    # If a HINV transcript has only one associated gene name/ENSG, then add those annotations to the transcript
    for hinv_transcript, annots in accum_annot.items():
        gene_names = filter(lambda x:x[0]=="RefSeq", annots)
        gencode_genes = filter(lambda x:x[0]=="GENCODE", annots)
        if (len(gencode_genes)==1 and len(gene_names)<=1):
            annot_per_transcript[hinv_transcript]["GENCODE"] = gencode_genes[0][1]
            if (len(gene_names)!=0):
                annot_per_transcript[hinv_transcript]["RefSeq"] = gene_names[0][1]

    return annot_per_transcript


def addAnnotaton(annot_per_transcript, unannot_gtf_gz, output_gtf_gz):
    reGTF = re.compile("^(chr.+\s+H-InvDB\s+mRNA\s+\d+\s+\d+\s+.+\s+.\s+.\s+.+)Name=(\S+\d);(.*)$")

    ip = gzip.open(unannot_gtf_gz, 'rb')
    op = gzip.open(output_gtf_gz, 'wb')

    pdb.set_trace()
    
    for line in ip:
        mo = reGTF.match(line)
        line = line.strip()
        if (mo != None):
            first_part, transcript_id, rest = mo.groups()

            if (annot_per_transcript[transcript_id].has_key("GENCODE")):
                rest += ";GENCODE=%s" % annot_per_transcript[transcript_id]["GENCODE"]
            if (annot_per_transcript[transcript_id].has_key("RefSeq")):
                rest += ";RefSeq=%s" % annot_per_transcript[transcript_id]["RefSeq"]

            op.write("%sName=%s;%s\n" % (first_part, transcript_id, rest))
        else:
            op.write("%s\n" % line)

    ip.close()
    op.close()


if (__name__ == "__main__"):
    intersect_file, unannot_gtf_gz, output_gtf_gz = sys.argv[1:]

    annot_per_transcript = compileAnnotation(intersect_file)
    addAnnotaton(annot_per_transcript, unannot_gtf_gz, output_gtf_gz)
    sys.exit(0)
