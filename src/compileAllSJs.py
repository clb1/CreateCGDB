#!/usr/bin/env python

import bz2
from operator import itemgetter
import sys
from parseExternalDatabase import readIsoformModels
from RNAIsoform import RNAIsoform


def getUniqueSJs(all_isoform_models):
    all_SJs = set()

    for chrom, isoforms in all_isoform_models.items():
        for isoform in isoforms:
            strand = isoform.getStrand()
            intron_tuples = isoform.getIntronTuples()
            for left_pos, right_pos in intron_tuples:
                all_SJs.add( (chrom, left_pos, right_pos, strand) )

    all_SJs = sorted(list(all_SJs), key=itemgetter(0,1,2,3))

    return all_SJs


def writeSJs(all_SJs, output_tsv):
    op = bz2.BZ2File(output_tsv, 'w')
    for tup in all_SJs:
        op.write("%s\t%d\t%d\t%s\n" % tup)
    op.close()


if (__name__ == "__main__"):
    isoform_model_databases, output_tsv = sys.argv[1:]
    
    all_isoform_models = readIsoformModels(isoform_model_databases, None)
    all_SJs = getUniqueSJs(all_isoform_models)
    writeSJs(all_SJs, output_tsv)

    sys.exit(0)
    
