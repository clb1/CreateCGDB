#!/usr/bin/env python

#import gzip
from collections import defaultdict
import cPickle
import numpy as np
from operator import itemgetter
import os
from string import find, maketrans, translate
import sys
from parseExternalDatabase import *
from RNAIsoformAnnotator import *
from RNAIsoform import RNAIsoform

#from Bio.Alphabet import IUPAC
#from Bio.Seq import Seq
import pysam

#import pybedtools
# Do this in main and put on command line
#pybedtools.set_bedtools_path("/raid/software/src/bedtools2/bin")

sys.path.append("/projects/CGDB2/other_software/liblinear-2.01/python")
sys.path.append("/projects/CGDB2/other_software/PolyAMotif")

import liblinearutil
from poly_a_class import PredictHMMClassifier


def readModels(models_directory, polyA_motifs):
    hmm_and_svm_models = {}
    for motif in polyA_motifs:
        hmms = cPickle.load(open("%s/%s_hmm.pkl" % (models_directory, motif), 'rb'))
        svm = liblinearutil.load_model("%s/%s.svm" %( models_directory, motif))
        hmm_and_svm_models[motif] = (hmms, svm)
    return hmm_and_svm_models


def orderChromosomes(all_isoform_models):
    ordered_chromosomes = []

    chromosomes = all_isoform_models.keys()

    numeric_chromosomes, alphnum_chromosomes = [], []
    for c in chromosomes:
        if (c[3:].isdigit()):
            numeric_chromosomes.append(int(c[3:]))
        else:
            alphnum_chromosomes.append(c[3:])
    numeric_chromosomes.sort()
    alphnum_chromosomes.sort()
    ordered_chromosomes = map(lambda x: "chr%s" % x, map(str, numeric_chromosomes) + alphnum_chromosomes)

    return ordered_chromosomes


def compilePolyASiteSequences(isoform_models, polyA_motifs): #, tempdir
    nucleotide_complement = maketrans("ACGT", "TGCA")
    number_models_w_polyA_motifs = 0
    distinct_region_labels = set()
    #import pdb
    #pdb.set_trace()
    
    polyA_regions_formatted_for_prediction = defaultdict(dict)
    for isoform in isoform_models:
        transcript_id = isoform.getTranscriptIDInSrcDB()
        chromosome = isoform.getChromosome()
        strand = isoform.getStrand()
        terminus_three_prime = isoform.getStrandwiseStop()

        # Perform sequence analysis and prediction sequence extraction on sequence oriented 5'->3'.
        nuc_seq = isoform.getThreePrimeTerminalRegionSequence(genome_fasta, 140)
        if (nuc_seq == None):
            continue
        else:
            nuc_seq = nuc_seq.upper()

        candidate_polyA_signal_region = nuc_seq[140-35:140-9]
        polyA_motifs_in_region = filter(lambda x: x in candidate_polyA_signal_region, polyA_motifs)

        if (len(polyA_motifs_in_region) > 0):
            number_models_w_polyA_motifs += 1
            for motif in polyA_motifs_in_region:
                indices = []
                ind = find(candidate_polyA_signal_region, motif)
                while(ind != -1):
                    indices.append( ind )
                    ind = find(candidate_polyA_signal_region, motif, ind+1)
                assert( len(indices) > 0 )

                for ind in indices:
                    # Cleavage site region is 10-30bp downstream of the end of the 6-bp polyA signal
                    chromosomal_cleavage_region_start = terminus_three_prime - 35 + ind + 5 + 10
                    chromosomal_cleavage_region_stop = terminus_three_prime - 35 + ind + 5 + 30

                    if (strand == '-'):
                        # Reflect start/stop positions about the isoform terminus position
                        new_start = terminus_three_prime - (chromosomal_cleavage_region_stop - terminus_three_prime)
                        new_stop = terminus_three_prime + (terminus_three_prime - chromosomal_cleavage_region_start)
                        chromosomal_cleavage_region_start = new_start
                        chromosomal_cleavage_region_stop = new_stop
                        #chromosomal_cleavage_region_start = terminus_three_prime + 35 - ind - 6 - 29 #10 + ind - 30
                        #chromosomal_cleavage_region_stop = terminus_three_prime + 35 -ind -6 -9 #10 + ind - 10

                    # Need 100bp to either side of polyA signal for prediction
                    seq_for_prediction = nuc_seq[140-35+ind-100:140-35+ind+5+100+1] # 100bp-motif-100bp

                    try:
                        assert (len(seq_for_prediction)==206)
                    except AssertionError:
                        import pdb
                        pdb.set_trace()
                        
                    region_label = (chromosome, chromosomal_cleavage_region_start, chromosomal_cleavage_region_stop, motif, strand)
                    if (not polyA_regions_formatted_for_prediction[motif].has_key(region_label)):
                        polyA_regions_formatted_for_prediction[motif][region_label] = seq_for_prediction
                        distinct_region_labels.add(region_label)
                        
    return (number_models_w_polyA_motifs, len(distinct_region_labels), polyA_regions_formatted_for_prediction)


def extractAndWriteGenomeRegions(pbt_merged, genome_fasta, op_fasta):
    for line in pbt_merged:
        chromosome, start, stop, strand = line.fields

        region_spec = "%s:%d-%s" % (chromosome, int(start)+1, stop)
        region_id = "%s:%s" % (region_spec, strand)
        nuc_seq_fasta = pysam.faidx(genome_fasta, region_spec)
        nuc_seq = ''.join(map(lambda x: x.strip(), nuc_seq_fasta[1:]))
        nuc_seq = nuc_seq.upper()
        nuc_seq = Seq(nuc_seq, IUPAC.unambiguous_dna)

        op_fasta.write(">%s\n%s\n" % (region_id, nuc_seq))
    

def predictAndWritePolyASites(polyA_regions_formatted_for_prediction, hmm_and_svm_models, output_bed):
    num_positive_predictions = 0
    d = 5
    k=20
    #C=1

    #import pdb
    #pdb.set_trace()
    tups_to_write = []
    for motif, D in polyA_regions_formatted_for_prediction.items():
        hmm_model_list, svm_model = hmm_and_svm_models[motif]
        test_data_names = []
        test_data = []
        for label, sequence in D.items():
            test_data_names.append(label)
            test_data.append( list(sequence) )
        test_lab = np.ones((len(test_data),), 'f')
                       
        print >> sys.stderr, "INFO: predicting for polyA signal motif %s (n=%d)..." % (motif, len(test_lab))
        test_predict, test_acc, predict_val = PredictHMMClassifier(test_data, test_lab, d, k, svm_model, hmm_model_list)
        sys.stderr.flush()
        sys.stdout.flush()
        
        # test_data_names is list of (chromosome, chromosomal_cleavage_region_left_boundary, chromosomal_cleavage_region_right_boundary, motif, strand)
        for tup, pred_val in zip(test_data_names, predict_val):
            if (pred_val > 0):
                tups_to_write.append( tup[0:4] + (1,) + tup[4:] )
                num_positive_predictions += 1
            else:
                tups_to_write.append( tup[0:4] + (0,) + tup[4:] )

    tups_to_write = sorted(tups_to_write, key=itemgetter(1))
    
    op = open(output_bed, 'w')
    for tup in tups_to_write:
        op.write("%s\t%d\t%d\t%s\t%d\t%s\n" % tup)
    op.close()

    return num_positive_predictions


if (__name__ == "__main__"):
    models_directory, polyA_motifs, genome_fasta, isoform_model_databases, output_bed = sys.argv[1:] # tempdir_root,
    polyA_motifs = polyA_motifs.split(',')
    
    tempdir = "%s/extractTerminiRegions_%s_%d" % (tempdir_root, os.getlogin(), os.getpid())
    os.mkdir(tempdir)
    pybedtools.set_tempdir(tempdir)
    hmm_and_svm_models = readModels(models_directory, polyA_motifs)
    
    annotator = RNAIsoformAnnotator()
    all_isoform_models = readIsoformModels(isoform_model_databases, annotator)

    ordered_chromosomes = orderChromosomes(all_isoform_models)

    for chromosome in ordered_chromosomes:
        print >> sys.stderr, "INFO: extracting 3' isoform termini regions on %s" % chromosome
        isoform_models = all_isoform_models[chromosome]

        num_models_w_polyA_motifs, num_distinct_polyA_locations, polyA_regions_formatted_for_prediction = compilePolyASiteSequences(isoform_models, polyA_motifs)
        print >> sys.stderr, "INFO: %d of %d isoforms had candidate polyA signal motifs, %d distinct." % (num_models_w_polyA_motifs, len(isoform_models), num_distinct_polyA_locations)
        print >> sys.stderr, "INFO: predicting actual polyA signal sequences"

        num_positive_predictions = predictAndWritePolyASites(polyA_regions_formatted_for_prediction, hmm_and_svm_models, output_bed)
        print >> sys.stderr, "INFO: predicted %d of %d polyA signals positively." % (num_positive_predictions, num_distinct_polyA_locations)
        
    os.rmdir(tempdir)

    sys.exit(0)
    
