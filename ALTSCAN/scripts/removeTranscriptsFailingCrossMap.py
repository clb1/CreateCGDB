#!/usr/bin/env python

import re
import sys
import pdb

if (__name__ == "__main__"):
    failed_loci = set()
    with open(sys.argv[1],'r') as ip:
        for line in ip:
            line = line.strip()
            if (line != ''):
                try:
                    locus_ID = line.rsplit('.',1)[0]
                except ValueError:
                    print >> sys.stderr, line
                failed_loci.add(locus_ID)
        
    unfiltered_IDs, filtered_IDs = set(), set()
    op_filtered = open("filtered_transcripts.gff", 'w')
    reTranscript = re.compile("transcript_id (\S+);")
    for line in sys.stdin:
        line = line.strip()
        mo = reTranscript.search(line)
        if (mo):
            transcript = mo.group(1)
            locus = transcript.rsplit('.',1)[0]
            if (locus not in failed_loci):
                unfiltered_IDs.add(transcript)
                print >> sys.stdout, line    
            else:
                filtered_IDs.add(transcript)
                op_filtered.write("%s\n" % line)
        else:
            print >> sys.stdout, line
    op_filtered.close()

    num_total_transcripts = len(filtered_IDs) + len(unfiltered_IDs)
    percent_filtered = 100.0 * float(len(filtered_IDs)) / float(num_total_transcripts)

    print >> sys.stderr, "INFO: filtered %d of %d total (%4.3f%%)" % (len(filtered_IDs), num_total_transcripts, percent_filtered)
    
    sys.exit(0)
    
