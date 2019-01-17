#!/usr/bin/env python

import sys
from parseExternalDatabase import *
from RNAIsoformAnnotator import *
from RNAIsoform import RNAIsoform


def formOutputBedName(db):
    output_bed = None
    if (".gff" in db):
        output_bed = db.split(".gff")[0] + ".bed"
    elif (".gtf" in db):
        output_bed = db.split(".gtf")[0] + ".bed"
    else:
        print >> sys.stderr, "ERROR: cannot convert isoform database name (%s) to corresponding BED file name" % db
        sys.exit(1)
        
    return output_bed


def writeBed12(all_isoform_models, output_bed):
    op = open(output_bed, 'w')

    isoform_models_for_chrom = all_isoform_models.values()
    for isoform_models in isoform_models_for_chrom:
        for isoform in isoform_models:
            isoform_fields = isoform.getAsBED12Fields("ID")
            op.write("%s\n" % "\t".join(isoform_fields))

    op.close()

    
if (__name__ == "__main__"):
    hgnc_data_tsv, isoform_model_databases = sys.argv[1:]

    annotator = RNAIsoformAnnotator(hgnc_data_tsv)
    
    for db in isoform_model_databases.split(','):
        output_bed = formOutputBedName(db)
        all_isoform_models = readIsoformModels(db, annotator)
        writeBed12(all_isoform_models, output_bed)    
        print >> sys.stderr, "INFO: converted %s to %s\n" % (db, output_bed)

    sys.exit(0)
    
