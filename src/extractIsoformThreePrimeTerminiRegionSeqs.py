#!/usr/bin/env python

import gzip
import os
import sys
from parseExternalDatabase import *
from RNAIsoformAnnotator import *
from RNAIsoform import RNAIsoform

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import pysam

import pybedtools
# Do this in main and put on command line
pybedtools.set_bedtools_path("/raid/software/src/bedtools2/bin")


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


def compileAndMergeRegionDefinitions(isoform_models, tempdir):
    unmerged_bed = "%s/first_exons_unsorted.bed" % tempdir
    op = open(unmerged_bed, "w")
    for isoform in isoform_models:
        chromosome = isoform.getChromosome()
        strand = isoform.getStrand()
        terminus_three_prime = isoform.getStrandwiseStop()
        bed_start = terminus_three_prime - 101
        bed_stop = terminus_three_prime + 100
        bed_line = "%s\t%d\t%d\tna\t0\t%s\n" % (chromosome, bed_start, bed_stop, strand)
        op.write(bed_line)
    op.close()
    
    pbt_unmerged = pybedtools.BedTool(unmerged_bed)
    pbt_unmerged_sorted = pbt_unmerged.sort()

    pbt_merged = pbt_unmerged_sorted.merge(s=True, c=6, o="distinct")
        
    os.remove(unmerged_bed)

    return pbt_merged


def extractAndWriteGenomeRegions(pbt_merged, genome_fasta, op_fasta):
    for line in pbt_merged:
        chromosome, start, stop, strand = line.fields

        region_spec = "%s:%d-%s" % (chromosome, int(start)+1, stop)
        region_id = "%s:%s" % (region_spec, strand)
        nuc_seq_fasta = pysam.faidx(genome_fasta, region_spec)
        nuc_seq = ''.join(map(lambda x: x.strip(), nuc_seq_fasta[1:]))
        nuc_seq = nuc_seq.upper()
        nuc_seq = Seq(nuc_seq, IUPAC.unambiguous_dna)
        if (strand == '-'):
            nuc_seq = nuc_seq.reverse_complement()

        op_fasta.write(">%s\n%s\n" % (region_id, nuc_seq))
    

if (__name__ == "__main__"):
    tempdir_root, genome_fasta, isoform_model_databases, output_fasta = sys.argv[1:]

    tempdir = "%s/extractTerminiRegions_%s_%d" % (tempdir_root, os.getlogin(), os.getpid())
    os.mkdir(tempdir)
    pybedtools.set_tempdir(tempdir)

    annotator = RNAIsoformAnnotator()
    all_isoform_models = readIsoformModels(isoform_model_databases, annotator)

    ordered_chromosomes = orderChromosomes(all_isoform_models)

    op_fasta = gzip.open(output_fasta, 'wb')
    for chromosome in ordered_chromosomes:
        print >> sys.stderr, "INFO: extracting 3' isoform termini regions on %s" % chromosome
        isoform_models = all_isoform_models[chromosome]
        pbt_merged = compileAndMergeRegionDefinitions(isoform_models, tempdir)
        extractAndWriteGenomeRegions(pbt_merged, genome_fasta, op_fasta)
    op_fasta.close()
    
    pybedtools.cleanup(remove_all=True)
    os.rmdir(tempdir)

    sys.exit(0)
    
