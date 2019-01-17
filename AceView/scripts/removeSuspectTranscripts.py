#!/usr/bin/env python

from operator import itemgetter
import re
import sys
import pdb

if (__name__ == "__main__"):
    mRNA_end_support = sys.argv[1]

    transcript_support = {}
    with open(mRNA_end_support, 'r') as ip:
        for line in ip:
            transcript_ID, mRNA_start_has_support, mRNA_end_has_support = line.strip().split("\t")
            transcript_support[transcript_ID] = (mRNA_start_has_support, mRNA_end_has_support)

    reTranscript = re.compile("^(chr.+)\s+AceView\s+exon\s+(\d+)\s+(\d+)\s+.\s+(\S)\s+\d\s+gene_id (.+); .+ transcript_id (.+); exon_number (\d+)")

    curr_transcript_lines = []
    curr_transcript_ID = None
    curr_chromosome = None
    curr_gene_id = None
    for line in sys.stdin:
        line = line.strip()

        mo = reTranscript.match(line)
        try:
            chromosome, exon_start, exon_stop, strand, gene_id, transcript_ID, exon_number = mo.groups()
        except AttributeError:
            print >> sys.stderr, line
            pdb.set_trace()

        if (curr_transcript_ID == None):
            curr_transcript_ID = transcript_ID
            curr_chromosome = chromosome
            curr_gene_id = gene_id

        if (transcript_ID == curr_transcript_ID):
            curr_transcript_lines.append( (line, int(exon_start), int(exon_stop), strand, int(exon_number)) )
            try:
                assert (chromosome == curr_chromosome), "Multiple chromsomes for transcript %s" % curr_transcript_ID
                assert (gene_id == curr_gene_id), "Multiple gene_id's for transcript %s" % curr_transcript_ID
            except AssertionError as ae:
                print >> sys.stderr, transcript_ID, curr_transcript_ID
                print >> sys.stderr, chromosome, curr_chromosome
                print >> sys.stderr, gene_id, curr_gene_id
                print >> sys.stderr, ae.message
                print >> sys.stderr, line
                pdb.set_trace()
        else:
            exon_starts = map(itemgetter(1), curr_transcript_lines)
            exon_stops = map(itemgetter(2), curr_transcript_lines)
            exon_strands = set(map(itemgetter(3), curr_transcript_lines))
            exon_numbers = map(itemgetter(4), curr_transcript_lines)

            #if (curr_transcript_ID.startswith("SLC25A11.cAug10")):
            #    print exon_numbers

            if (exon_numbers[0] != 1):
                print >> sys.stderr, "%s skipped, first exon not exon_number 1" % curr_transcript_ID
            elif (len(exon_strands) > 1):
                print >> sys.stderr, "\n".join(map(itemgetter(0), curr_transcript_lines))
                print >> sys.stderr, "%s skipped, exons from both strands" % curr_transcript_ID
            else:
                if (all(map(lambda (a,b):a-b==1, zip(exon_numbers[1:], exon_numbers[:-1])))):
                    assert (len(exon_starts) == len(exon_stops))
                    mRNA_start = min(exon_starts)
                    mRNA_stop = max(exon_stops)
                    curr_strand = list(exon_strands)[0]
                    if (transcript_support.has_key(curr_transcript_ID)):
                        support_5p, support_3p = transcript_support[curr_transcript_ID]
                    else:
                        support_5p, support_3p = ('U','U')
                    tup = (curr_chromosome, mRNA_start, mRNA_stop, curr_strand, curr_gene_id, curr_transcript_ID, support_5p, support_3p, len(exon_starts))
                    print >> sys.stdout, "%s\tAceView\tmRNA\t%d\t%d\t.\t%s\tgene_id %s; transcript_id %s; mRNA_start_support %s; mRNA_end_support %s; number_exons %d" % tup
                    for l in map(itemgetter(0), curr_transcript_lines):
                        print >> sys.stdout, l
                else:
                    print >> sys.stderr, "%s skipped, non-contiguous exon numbering" % curr_transcript_ID

            curr_transcript_lines = [ (line, int(exon_start), int(exon_stop), strand, int(exon_number)) ]
            curr_transcript_ID = transcript_ID
            curr_chromosome = chromosome
            curr_gene_id = gene_id


    sys.exit(0)
    
