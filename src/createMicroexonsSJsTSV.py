#!/usr/bin/env python

from collections import defaultdict
from operator import itemgetter
import os
import sys
from parseExternalDatabase import readIsoformModels
from RNAIsoform import RNAIsoform

import pybedtools


def getUniqueTranscripts(all_isoform_models):
    all_transcripts = defaultdict(list)

    for chrom, isoforms in all_isoform_models.items():
        for isoform in isoforms:
            isoform.setLengthAndExonTermini()
            transcript_id = isoform.getTranscriptIDInSrcDB()
            strand = isoform.getStrand()
            left_pos, right_pos = isoform.getStartStop()
            tup = (chrom, left_pos, right_pos, strand)
            all_transcripts[tup].append(transcript_id)

    return all_transcripts


def writeTranscriptBoundaries(all_transcripts):
    bed_tuples = []
    for tup, transcript_ids in all_transcripts.items():
        for transcript_id in transcript_ids:
            bed_tuples.append( (tup[0], tup[1], tup[2], transcript_id, tup[3]) )
        
    bed_tuples = sorted(bed_tuples, key=itemgetter(0,1))

    boundaries_bed = "temp_transcript_boundaries.bed"
    op = open(boundaries_bed, 'w')
    for tup in bed_tuples:
        op.write("%s\t%d\t%d\t%s\t1\t%s\n" % tup)
    op.close()

    return boundaries_bed


def identifyParentTranscripts(microexons_bed, boundaries_bed, all_transcripts):
    """Returns all coordinates as 1-based."""
    exclude_existing_microexons = True
    
    microexons_per_transcript = defaultdict(set)
    known_exons = set()
    
    microexons = pybedtools.BedTool(microexons_bed)
    transcripts = pybedtools.BedTool(boundaries_bed)

    b = transcripts.intersect(microexons, loj=True, stream=True)
    for line in b:
        fields = line.fields
        if (fields[6] != '.'):
            me_name = fields[9]
            microexons_per_transcript[fields[3]].add( (fields[6], int(fields[7])+1, int(fields[8]), fields[11]) )
            known_exons.add( (fields[6], int(fields[7])+1, int(fields[8])) )

    # Remove microexons that already exist in any isoform, as they are already "accounted for".
    if (exclude_existing_microexons):
        exon_tuples_to_remove = set()
        for isoforms in all_isoform_models.values():
            for isoform in isoforms:
                transcript_id = isoform.getTranscriptIDInSrcDB()
                if (microexons_per_transcript.has_key(transcript_id)):
                    exon_tuples = isoform.getExonsAsTuples()
                    for exon_tuple in exon_tuples:
                        #if (exon_tuple == ("chr8", 28429378, 28429403)):
                        #    print >> sys.stderr, transcript_id
                        if (exon_tuple in known_exons):
                            exon_tuples_to_remove.add( exon_tuple )

        transcript_ids_to_remove = set()
        for transcript_id, exon_tuples in microexons_per_transcript.items():
            ets_to_remove = set()
            for et in exon_tuples:
                trunc_exon_tuple = (et[0], et[1], et[2])
                if (trunc_exon_tuple in exon_tuples_to_remove):
                    ets_to_remove.add(et)

            if (len(ets_to_remove) > 0):
                microexons_per_transcript[transcript_id] = exon_tuples - ets_to_remove
                if (len(microexons_per_transcript[transcript_id]) == 0):
                    transcript_ids_to_remove.add(transcript_id)

        for transcript_id in transcript_ids_to_remove:
            del microexons_per_transcript[transcript_id]

    return microexons_per_transcript

    
def createNewSJs(microexons_per_transcript, all_isoform_models):
    new_SJs = set()

    for chrom, isoforms in all_isoform_models.items():
        for isoform in isoforms:
            transcript_id = isoform.getTranscriptIDInSrcDB()
            if (microexons_per_transcript.has_key(transcript_id) and isoform.getNumberOfExons() > 1):
                transcript_strand = isoform.getStrand()
                exon_5p = isoform.getExonFivePrimeTerminiPositions()
                exon_3p = isoform.getExonThreePrimeTerminiPositions()

                for me_chrom, me_start, me_stop, me_strand in microexons_per_transcript[transcript_id]:
                    assert (chrom == me_chrom)
                    if (transcript_strand == '-' and me_strand in ['-','.']):
                        for pos_3p in exon_3p:
                            if (pos_3p > me_stop):
                                new_SJs.add( (chrom, me_stop+1, pos_3p-1, transcript_strand) )
                        for pos_5p in exon_5p:
                            if (pos_5p < me_start):
                                new_SJs.add( (chrom, pos_5p+1, me_start-1, transcript_strand) )
                    elif (transcript_strand == '+' and me_strand in ['+','.']):
                        for pos_3p in exon_3p:
                            if (pos_3p < me_start):
                                new_SJs.add( (chrom, pos_3p+1, me_start-1, transcript_strand) )
                        for pos_5p in exon_5p:
                            if (pos_5p > me_stop):
                                new_SJs.add( (chrom, me_stop+1, pos_5p-1, transcript_strand) )

    return new_SJs


def writeTSV(new_SJs, output_tsv):
    new_SJs = sorted(list(new_SJs), key=itemgetter(0, 1))
    with open(output_tsv, 'w') as op:
        for tup in new_SJs:
            op.write("%s\t%d\t%d\t%s\n" % tup)


if (__name__ == "__main__"):
    microexons_bed, isoform_model_databases, output_tsv = sys.argv[1:]
    
    all_isoform_models = readIsoformModels(isoform_model_databases, None)
    all_transcripts = getUniqueTranscripts(all_isoform_models)
    boundaries_bed = writeTranscriptBoundaries(all_transcripts)

    microexons_per_transcript = identifyParentTranscripts(microexons_bed, boundaries_bed, all_transcripts)

    new_SJs = createNewSJs(microexons_per_transcript, all_isoform_models)

    writeTSV(new_SJs, output_tsv)

    # Cleanup
    os.unlink(boundaries_bed)
    
    sys.exit(0)
    
