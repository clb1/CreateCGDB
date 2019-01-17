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
        
        aceview_exon_start = fields[12]
        aceview_exon_stop = fields[13]
        aceview_annot = fields[17]
        
        if (gencode_exon_start==aceview_exon_start or gencode_exon_stop==aceview_exon_stop):
            gencode_annot_elems = gencode_annot.split()
            gencode_annot_dict = dict([(gencode_annot_elems[i],gencode_annot_elems[i+1]) for i in xrange(0,len(gencode_annot_elems),2)])
            aceview_annot_elems = aceview_annot.split()
            aceview_annot_dict = dict([(aceview_annot_elems[i],aceview_annot_elems[i+1]) for i in xrange(0,len(aceview_annot_elems),2)])

            aceview_gene = aceview_annot_dict['gene_id'][0:-1]
            aceview_transcript = aceview_annot_dict['transcript_id'][0:-1]
            gencode_gene = gencode_annot_dict['gene_id'][1:-2] # Has " characters
            accum_annot[aceview_transcript].add( ("GENCODE", gencode_gene) )
            if (gencode_annot_dict.has_key("gene_name")):
                accum_annot[aceview_transcript].add( ("RefSeq", gencode_annot_dict["gene_name"][1:-2]) )
            if (not aceview_gene[0].islower() and 'and' not in aceview_gene):
                accum_annot[aceview_transcript].add( ("AceView", aceview_gene) )
    ip.close()

    # If a AceView transcript has only one associated gene name/ENSG, then add those annotations to the transcript
    for aceview_transcript, annots in accum_annot.items():
        gene_ids = filter(lambda x:x[0]=="AceView", annots)
        try:
            assert (len(gene_ids) <= 1), "Multiple Aceview gene_id's for same transcript"
        except AssertionError as ae:
            pdb.set_trace()
        gene_id = gene_ids[0][1] if (len(gene_ids)==1) else None

        gene_names = filter(lambda x:x[0]=="RefSeq", annots)
        gencode_genes = filter(lambda x:x[0]=="GENCODE", annots)
        if (len(gencode_genes)==1 and len(gene_names)<=1):
            annot_per_transcript[aceview_transcript]["GENCODE"] = gencode_genes[0][1]
            if (len(gene_names)!=0):
                gene_name = gene_names[0][1]
                #if (gene_id != None and gene_name != gene_id):
                #    pdb.set_trace()
                annot_per_transcript[aceview_transcript]["RefSeq"] = gene_name

    return annot_per_transcript


def addAnnotaton(annot_per_transcript, unannot_gtf_gz, output_gtf_gz):
    reGTF = re.compile("^(chr.+\s+AceView\s+mRNA\s+\d+\s+\d+\s+.+\s+.\s+gene_id \S+; transcript_id) (\S+);(.+)$")

    ip = gzip.open(unannot_gtf_gz, 'rb')
    op = gzip.open(output_gtf_gz, 'wb')

    for line in ip:
        mo = reGTF.match(line)
        line = line.strip()
        if (mo != None):
            first_part, transcript_id, rest = mo.groups()

            if (annot_per_transcript[transcript_id].has_key("GENCODE")):
                rest += "; GENCODE %s" % annot_per_transcript[transcript_id]["GENCODE"]
            if (annot_per_transcript[transcript_id].has_key("RefSeq")):
                rest += "; RefSeq %s" % annot_per_transcript[transcript_id]["RefSeq"]

            op.write("%s %s;%s\n" % (first_part, transcript_id, rest))
        else:
            op.write("%s\n" % line)

    ip.close()
    op.close()


if (__name__ == "__main__"):
    intersect_file, unannot_gtf_gz, output_gtf_gz = sys.argv[1:]

    annot_per_transcript = compileAnnotation(intersect_file)
    addAnnotaton(annot_per_transcript, unannot_gtf_gz, output_gtf_gz)
    sys.exit(0)
