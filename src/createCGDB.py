#!/usr/bin/env python
from __future__ import print_function

from collections import defaultdict, OrderedDict
from itertools import chain, combinations, product
from operator import itemgetter, attrgetter, xor, methodcaller
import cPickle
from pyfaidx import Fasta
import gzip
import networkx as nx
import os
import shutil
import sys
from parseExternalDatabase import *
from RNAIsoformAnnotator import *
from RNAIsoform_python2 import Exon, RNAIsoform
from IsoformGroup import IsoformGroup
from TSS import TSS

import pdb

import pybedtools

def orderChromosomes(all_TSSs, all_isoform_models, filter_nonstandard_contigs=True):
    ordered_chromosomes = []

    if (filter_nonstandard_contigs):
        for contig_name in all_TSSs.keys():
            if (contig_name.endswith("_alt") or contig_name.endswith("_random") or contig_name.startswith("chrUn_")):
                del all_TSSs[contig_name]

        for contig_name in all_isoform_models.keys():
            if (contig_name.endswith("_alt") or contig_name.endswith("_random") or contig_name.startswith("chrUn_")):
                del all_isoform_models[contig_name]

    assert (set(all_TSSs.keys()) <= set(all_isoform_models.keys()))

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


def standardizePolyASites(pAs_and_exons):
    """'pAs_and_exons' is list of PolyA objects and last exon Exon objects.
    After this function, all last exons will have been adjusted to standardized
    polyA sites for CGDB."""

    # Normal range for a polyA site downstream of a polyA motif is 10-30bp.
    # Added 3bp of leeway.
    pAs_range_limit = 23
    
    pAs_and_exon_3ps = []
    pAs_and_exons_strands = set()

    for c in filter(lambda x: isinstance(x, Exon), pAs_and_exons):
        if (c.isThreePrimeComplete()):
            exon_3p_pos = c.getThreePrimePosition()
            pAs_and_exons_strands.add(c.getStrand())
            pAs_and_exon_3ps.append( (exon_3p_pos, c) )

    for c in filter(lambda x: isinstance(x, PolyA), pAs_and_exons):
        site_pos = c.getSitePosition()
        pAs_and_exons_strands.add(c.getStrand())
        pAs_and_exon_3ps.append( (site_pos, c) )
        
    assert (len(pAs_and_exons_strands)==1)
        
    strand = list(pAs_and_exons_strands)[0]

    # Logic below processes polyA sites in a 3'->5' manner, regardless of strand
    remaining_pAs_and_exon_3ps_to_group = sorted(pAs_and_exon_3ps, key=itemgetter(0), reverse=(strand=='+'))

    # Greedy approach.
    # Accumulate groups of distinct sites that all could be feasibly associated with polyA signal (ie ~10-30bp downstream).
    while (len(remaining_pAs_and_exon_3ps_to_group) > 0):
        all_site_groupings = []
        for i, init_tup in enumerate(remaining_pAs_and_exon_3ps_to_group):
            curr_site_group = set([init_tup[0]])
            for tup in remaining_pAs_and_exon_3ps_to_group[i+1:]:
                if (abs(tup[0] - init_tup[0]) + 1 <= pAs_range_limit):
                    curr_site_group.add(tup[0])
            all_site_groupings.append( curr_site_group )

        # Determine the polyA site grouping with the largest number of members
        index_of_largest_group = 0
        for k in xrange(1, len(all_site_groupings)):
            if (len(all_site_groupings[k]) > len(all_site_groupings[index_of_largest_group])):
                index_of_largest_group = k

        # Standardize polyA site for exons in largest group.
        # Remove all PolyA and Exon objects from list of remaining cases.
        sites_of_largest_group = all_site_groupings[index_of_largest_group]
        standardized_polyA_pos = max(sites_of_largest_group) if (strand=='+') else min(sites_of_largest_group)
        cases_still_remaining = []
        for site_pos, c in remaining_pAs_and_exon_3ps_to_group:
            if (site_pos in sites_of_largest_group):
                if (isinstance(c, Exon)):
                    c.resetStrandwiseStop(standardized_polyA_pos)
            else:
                cases_still_remaining.append( (site_pos, c) )

        remaining_pAs_and_exon_3ps_to_group = cases_still_remaining


def standardizePromoterRegion2(TSSs_and_exons, genome_seq):
    """'TSSs_and_exons' is list of TSSs and first exons. After this function,
    all first exons for the promoter will have been adjusted to standardized
    transcription start sites for CGDB.

    Precondition: 'TSSs_and_exons' is a list longer than 1 and it contains at least
    one exon.
    """
    max_allowed_new_orf_len = 24 # Length is the approximately the number of amino acids in a new signal peptide plus the minimum antibody epitope.
    
    exons_w_adjusted_starts = set()
    all_fantom = {"robust":[], "permissive":[]}
    candidate_standard_start_positions = set()

    # Collect candidate standard start positions
    TSSs_and_exons_strands = set()
    for e in filter(lambda x: isinstance(x, Exon), TSSs_and_exons):
        start_pos = e.getFivePrimePosition()
        if (e.isFivePrimeComplete() != False):
            candidate_standard_start_positions.add(start_pos)
        TSSs_and_exons_strands.add(e.getStrand())

    for t in filter(lambda x: isinstance(x, TSS), TSSs_and_exons):
        TSSs_and_exons_strands.add(t.getStrand())
        assert (t.isSrcLabel("FANTOM", "robust")), "Should only be using robust FANTOM TSSs"
        tss_five_prime, tss_three_prime = t.getFivePrimeThreePrimePositions()
        tss_mid = (tss_five_prime + tss_three_prime)/2
        candidate_standard_start_positions.update([tss_five_prime, tss_mid, tss_three_prime])
            
    assert (len(TSSs_and_exons_strands)==1)
    strand = list(TSSs_and_exons_strands)[0]

    # Associate with each candidate start position those exons with which it is compatible
    G = nx.Graph()
    for start_pos, exon in product(candidate_standard_start_positions, filter(lambda x: isinstance(x, Exon), TSSs_and_exons)):
        exon_5p = exon.getStrandwiseStart()
        exon_3p = exon.getStrandwiseStop()

        # Prevent truncation of already incomplete 5' ends
        if (exon.isFivePrimeComplete() == False and exon.doesOverlapPosition(start_pos)):
            continue

        if (abs(start_pos-exon_5p)<=3):
            if ((strand=='+' and (start_pos<exon_3p or exon_5p==exon_3p)) or (strand=='-' and (start_pos>exon_3p or exon_5p==exon_3p))):
                G.add_edge(start_pos, exon)
        else:
            start_is_upstream = (strand=='+' and start_pos<=exon_5p) or (strand=='-' and start_pos>=exon_5p)
            new_start_modifies_orf = True
            new_start_severely_truncates = True
            if (start_is_upstream):
                new_start_modifies_orf = exon.fivePrimeExtensionCreatesORF(start_pos, genome_seq, max_allowed_new_orf_len)
                new_start_severely_truncates = False

            elif (exon.doesOverlapPosition(start_pos)):
                new_start_modifies_orf = exon.fivePrimeContractionRemovesORF(start_pos, genome_seq)
                new_start_severely_truncates = float(abs(start_pos-exon_3p))/float(abs(exon_5p-exon_3p)) < 0.5

            if (not (new_start_modifies_orf or new_start_severely_truncates)):
                G.add_edge(start_pos, exon)

    # Greedily standardize exon starts for those exons that can be standardized
    while (G.size() > 0):
        start_pos_and_exons = []
        for start_pos in filter(lambda x: isinstance(x, int), G.nodes()):
            if (G.degree(start_pos) > 0):
                # Accumulate the exon length modification amounts if this start_pos used
                sum_abs_diff = 0
                unique_exons = set()
                for exon in G.neighbors(start_pos):
                    exon_5p = exon.getStrandwiseStart()
                    exon_3p = exon.getStrandwiseStop()
                    start_is_upstream = (strand=='+' and start_pos<=exon_5p) or (strand=='-' and start_pos>=exon_5p)
                    if (e.isFivePrimeComplete() != False or start_is_upstream):
                        unique_exons.add( (exon_5p,exon_3p) )
                        sum_abs_diff += abs(exon_5p-start_pos)
                if (len(unique_exons) > 0):
                    start_pos_and_exons.append( (start_pos, (sum_abs_diff+1.0)/float(len(unique_exons)), sum_abs_diff, len(unique_exons)) )
        start_pos_and_exons.sort(key=itemgetter(1))  # Convervatively groups TSSs

        new_start_pos = start_pos_and_exons[0][0]

        for exon in G.neighbors(new_start_pos):
            assert (isinstance(exon, Exon))
            exon.resetStrandwiseStart(new_start_pos)
            G.remove_node(exon)
        G.remove_node(new_start_pos)


def standardizePromoterRegion(TSSs_and_exons, genome_seq): # closest_upstream_tss, 
    """'TSSs_and_exons' is list of TSSs and first exons. After this function,
    all first exons for the promoter will have been adjusted to standardized
    transcription start sites for CGDB.

    Precondition: 'TSSs_and_exons' is a list longer than 1 and it contains at least
    one exon.
    """
    max_allowed_new_orf_len = 8 # Length is the number of amino acids. Is just larger than the minimum antibody epitope and min length of a primer.
    
    exons_w_adjusted_starts = set()
    all_fantom = {"robust":[], "permissive":[]}

    TSSs_and_exons_at_start_positions = defaultdict(list)
    TSSs_and_exons_strands = set()
    for c in filter(lambda x: isinstance(x, Exon), TSSs_and_exons):
        start_pos = c.getFivePrimePosition()
        TSSs_and_exons_strands.add(c.getStrand())
        TSSs_and_exons_at_start_positions[start_pos].append(c)

    # Add a TSS to any exon start position that it overlaps
    all_exon_starts = TSSs_and_exons_at_start_positions.keys()
    for c in filter(lambda x: isinstance(x, TSS), TSSs_and_exons):
        if (c.isSrcLabel("FANTOM", "robust")):
            all_fantom["robust"].append(c)
        else:
            all_fantom["permissive"].append(c)
        does_overlap_an_exon_start = False
        for exon_start_pos in all_exon_starts:
            if (c.doesOverlap(exon_start_pos)):
                does_overlap_an_exon_start = True
                TSSs_and_exons_at_start_positions[exon_start_pos].append(c)        
        if (not does_overlap_an_exon_start):
            start_pos = c.getFivePrimePosition()
            TSSs_and_exons_at_start_positions[start_pos].append(c)
            
    assert (len(TSSs_and_exons_strands)==1)
    strand = list(TSSs_and_exons_strands)[0]
    
    sorted_TSSs_and_exons = TSSs_and_exons_at_start_positions.items()
    sorted_TSSs_and_exons = sorted(sorted_TSSs_and_exons, key=lambda x:x[0], reverse=(strand=='-'))

    # Strategy: Find the most upstream TSS or other exon start to which each exon can be extended,
    #           unless it is already at a robust TSS or is already 5' complete.
    for i, (pos,TSSs_and_exons_list) in enumerate(sorted_TSSs_and_exons):
        exons = filter(lambda x: isinstance(x, Exon), TSSs_and_exons_list)

        if (len(exons) > 0):
            for fantom5_label in ["robust"]: # , "permissive"
                fantom_tss_at_pos = filter(lambda y: isinstance(y, TSS) and y.isSrcLabel("FANTOM", fantom5_label), TSSs_and_exons_list)

                # Try to standardize starts at this position based on a FANTOM TSS that overlaps the exon start position
                if (len(fantom_tss_at_pos) != 0):
                    assert (len(fantom_tss_at_pos) == 1)
                        
                    tss_five_prime, tss_three_prime = fantom_tss_at_pos[0].getFivePrimeThreePrimePositions()
                    new_tss = pos if (any(map(lambda exon: exon.isFivePrimeComplete(), exons))) else (tss_five_prime + tss_three_prime)/2
                    for e in exons:
                        if (e.isFivePrimeComplete() != False):
                            e.resetStrandwiseStart(new_tss)
                            exons_w_adjusted_starts.add(e)
                    break

                # If couldn't standardize based on overlapping FANTOM TSS, try to standardize starts at this position based
                # on a FANTOM TSS that is within its halfwidth (either upstream or downstream) of the exon start position
                #if (len(all_fantom[fantom5_label]) > 0):
                #    abs_dist_and_tss_tuples = filter(lambda y: y != None, map(lambda x: x.isWithinHalfwidth(pos), all_fantom[fantom5_label]))
                #    if (len(abs_dist_and_tss_tuples) > 0):
                #        abs_dist, tss_relative_direction, closest_tss = sorted(abs_dist_and_tss_tuples, key=lambda w:abs(w[0]))[0]
                #        tss_five_prime, tss_three_prime = closest_tss.getFivePrimeThreePrimePositions()
                #        for e in exons:
                #            if (e.isFivePrimeComplete()): # tss_relative_direction != "downstream" or 
                #                e.resetStrandwiseStart(tss_five_prime)
                #                delete e.setIfMostDownstreamVariableStart(tss_three_prime, fantom5_label)
                #                delete e.setIfMostDownstreamVariableStart(pos, fantom5_label)
                #                exons_w_adjusted_starts.add(e)
                #        break

                # Try to standardize based on an upstream exon start or FANTOM TSS
                for j in xrange(i):
                    # Don't continue if there are no exons for which to attempt extensions
                    if (not any(map(lambda x: x.isFivePrimeComplete() and x not in exons_w_adjusted_starts, exons))):
                        break
                    
                    upstream_fantom_tss = filter(lambda y: isinstance(y, TSS) and y.isSrcLabel("FANTOM", fantom5_label), sorted_TSSs_and_exons[j][1])
                    upstream_exons = filter(lambda y: isinstance(y, Exon) and not y in exons_w_adjusted_starts, sorted_TSSs_and_exons[j][1])
                    
                    if (len(upstream_fantom_tss) > 0):
                        assert (len(upstream_fantom_tss) == 1)
                        tss_five_prime, tss_three_prime = upstream_fantom_tss[0].getFivePrimeThreePrimePositions()

                    if (len(upstream_exons) > 0):
                        upstream_start = sorted_TSSs_and_exons[j][0]    

                    for e in exons:
                        if (e.isFivePrimeComplete()==False and e not in exons_w_adjusted_starts):
                            # Try to standardize exon at this position based on an upstream first exon
                            if (len(upstream_exons) > 0):
                                extension_creates_orf = e.fivePrimeExtensionCreatesORF(upstream_start, genome_seq, max_allowed_new_orf_len)
                                if (not extension_creates_orf):
                                    e.resetStrandwiseStart(upstream_start)
                                    exons_w_adjusted_starts.add(e)

                            # ...or on an upstream FANTOM TSS
                            elif (len(upstream_fantom_tss) > 0):
                                extension_creates_orf = e.fivePrimeExtensionCreatesORF(tss_three_prime, genome_seq, max_allowed_new_orf_len)
                                if (not extension_creates_orf):
                                    e.resetStrandwiseStart(tss_three_prime)
                                    exons_w_adjusted_starts.add(e)

    
def groupPolyASites(isoform_models, polyAsites, tempdir):
    lookup = {}
    G = nx.Graph()
    
    # Want:
    #       1) pairs of last exons that overlap
    #            - Use 'intersect' on last_exons.bed vs last_exons.bed
    #       2) last exon stops that are within polyA site variance range of a polyA site
    #            - Use 'slop' with 'b=0.5, pct=True' to expand a polyAsite.bed file, then intersect it with an exon_stops.bed
    #       3?) Nearest downstream polyA site of each last exon stop
    #            - Used 'closest' with 'ignore overlaps' and 'ignore upstream' options set.
    #            - Results must be used carefully so that a polyA site that is much too far downstream isn't later used.

    last_exons_bed = "%s/last_exons_unsorted.bed" % tempdir    
    last_exon_stops_bed = "%s/last_exon_stops_unsorted.bed" % tempdir
    pAs_bed = "%s/polyAsite_unsorted.bed" % tempdir

    op_exons_bed = open(last_exons_bed, 'w')
    for isoform_model in isoform_models:
        (ID, bed_line) = isoform_model.getLastExonForLabeledBED()
        lookup[ID] = isoform_model.getLastExon()
        op_exons_bed.write("%s\n" % bed_line)
    op_exons_bed.close()
    
    op_exon_stops_bed = open(last_exon_stops_bed, 'w')
    for isoform_model in isoform_models:
        (ID, bed_line) = isoform_model.getLastExonStrandwiseStopForLabeledBED()
        op_exon_stops_bed.write("%s\n" % bed_line)
    op_exon_stops_bed.close()

    op_pAs_bed = open(pAs_bed, 'w')
    for pAs in polyAsites:
        (ID, bed_line) = pAs.getPolyASiteForLabeledBED(slop=25)
        lookup[ID] = pAs
        op_pAs_bed.write("%s\n" % bed_line)
    op_pAs_bed.close()

    # Sort the BED files
    pbt_last_exons = pybedtools.BedTool(last_exons_bed)
    pbt_last_exons_sorted = pbt_last_exons.sort()
    
    pbt_last_exon_stops = pybedtools.BedTool(last_exon_stops_bed)
    pbt_last_exon_stops_sorted = pbt_last_exon_stops.sort()

    pbt_pAs = pybedtools.BedTool(pAs_bed)
    pbt_pAs_sorted = pbt_pAs.sort()

    # Connect overlapping last exons
    a = pbt_last_exons_sorted.intersect(pbt_last_exons_sorted, s=True, sorted=True, wo=True, stream=True)
    for line in a:
        if (line.fields[3] != line.fields[9]):
            elem_1 = lookup[ line.fields[3] ]
            elem_2 = lookup[ line.fields[9] ]
            G.add_edge(elem_1, elem_2)

    # Connect last exons with overlapping and nearby (downstream) polyAsites
    b = pbt_pAs_sorted.intersect(pbt_last_exon_stops_sorted, s=True, wo=True, stream=True)
    for line in b:
        elem_1 = lookup[ line.fields[3] ]
        elem_2 = lookup[ line.fields[9] ]
        if (not G.has_edge(elem_1, elem_2)):
            G.add_edge(elem_1, elem_2)

    print("INFO: standardizing polyA sites", file=sys.stderr)
    for pAs_and_exons in nx.connected_components(G):
        if (len(pAs_and_exons)>1 and any(map(lambda c: isinstance(c, Exon) and c.isThreePrimeComplete(), pAs_and_exons))):
            standardizePolyASites(pAs_and_exons)

    # Cleanup
    G.clear()
    G = None
    #os.remove(last_exons_bed)
    #os.remove(last_exon_stops_bed)
    #os.remove(pAs_bed)
    

def groupTranscriptionStarts(isoform_models, TSSs, genome_seq, tempdir):
    lookup = {}
    G = nx.Graph()

    print("INFO: grouping transcription start sites", file=sys.stderr)
    
    # Want:
    #       1) pairs of 1st exons that overlap
    #            - Use 'intersect' on first_exons.bed vs first_exons.bed
    #       2) 1st exon starts that are within halfwidth of TSS
    #            - Use 'slop' with 'b=0.5, pct=True' to expand a TSS.bed file, then intersect it with an exon_starts.bed
    #       3?) Nearest upstream TSS of each exon start
    #            - Used 'closest' with 'ignore overlaps' and 'ignore downstream' options set.
    #            - Results must be used carefully so that a TSS that is much too far upstream isn't later used.
    
    first_exons_bed = "%s/first_exons_unsorted.bed" % tempdir    
    first_exon_starts_bed = "%s/first_exon_starts_unsorted.bed" % tempdir
    tss_bed = "%s/tss_unsorted.bed" % tempdir

    op_exons_bed = open(first_exons_bed, 'w')
    for isoform_model in isoform_models:
        (ID, bed_line) = isoform_model.getFirstExonForLabeledBED()
        lookup[ID] = isoform_model.getFirstExon()
        op_exons_bed.write("%s\n" % bed_line)
    op_exons_bed.close()

    op_exon_starts_bed = open(first_exon_starts_bed, 'w')
    for isoform_model in isoform_models:
        (ID, bed_line) = isoform_model.getFirstExonStrandwiseStartForLabeledBED()
        op_exon_starts_bed.write("%s\n" % bed_line)
    op_exon_starts_bed.close()

    op_tss_bed = open(tss_bed, 'w')
    for tss in TSSs:
        (ID, bed_line) = tss.getTSSForLabeledBED(slop=0.5)
        lookup[ID] = tss
        op_tss_bed.write("%s\n" % bed_line)
    op_tss_bed.close()

    # Sort the BED files
    pbt_first_exons = pybedtools.BedTool(first_exons_bed)
    pbt_first_exons_sorted = pbt_first_exons.sort()
    
    pbt_first_exon_starts = pybedtools.BedTool(first_exon_starts_bed)
    pbt_first_exon_starts_sorted = pbt_first_exon_starts.sort()

    pbt_tss = pybedtools.BedTool(tss_bed)
    pbt_tss_sorted = pbt_tss.sort()

    # Connect overlapping first exons
    for line in pbt_first_exons_sorted.intersect(pbt_first_exons_sorted, s=True, sorted=True, wo=True, stream=True):
        if (line.fields[3] != line.fields[9]):
            elem_1 = lookup[ line.fields[3] ]
            elem_2 = lookup[ line.fields[9] ]
            G.add_edge(elem_1, elem_2)

    # Connect first exons and TSS's (for exons within halfwidth distance of the TSS)
    for line in pbt_tss_sorted.intersect(pbt_first_exon_starts_sorted, s=True, wo=True, stream=True):
        elem_1 = lookup[ line.fields[3] ]
        elem_2 = lookup[ line.fields[9] ]
        if (not G.has_edge(elem_1, elem_2)):
            G.add_edge(elem_1, elem_2)

    # Connect exons with nearest upstream TSS that is not overlapping
    #c = pbt_first_exon_starts_sorted.closest(pbt_tss_sorted, N=True, D="a", s=True, id=True, io=True, stream=True)
    #closest_upstream_tss = {}
    #for line in c:
    #    if (line.fields[6] == '.'):
    #        print("Skipped c: %s" % line
    #    else:
    #        elem_a = line.fields[3]
    #        elem_b = line.fields[9]
    #        closest_upstream_tss[lookup[elem_a]] = lookup[elem_b]
            
    #pybedtools.cleanup(remove_all=True)
    #print pbt_first_exons_sorted.fn
    #print pbt_first_exon_starts_sorted.fn
    #print pbt_tss_sorted.fn
    
    print("INFO: standardizing promoter regions", file=sys.stderr)
    for TSSs_and_exons in nx.connected_components(G):
        if (len(TSSs_and_exons)>1 and any(map(lambda c: isinstance(c, Exon), TSSs_and_exons))):
            standardizePromoterRegion2(TSSs_and_exons, genome_seq) # closest_upstream_tss, 

    # Cleanup
    G.clear()
    G = None
    #os.remove(first_exons_bed)
    #os.remove(first_exon_starts_bed)
    #os.remove(tss_bed)


def reportIsoformTerminiLengthAdjustments(isoform_models):
    isoform_length_adjustments = []
    for isoform in isoform_models:
        first_exon_length_adjustment_5p, last_exon_length_adjustment_3p = isoform.getTerminiLengthAdjustments()
        total_abs_length_adj = abs(first_exon_length_adjustment_5p) + abs(last_exon_length_adjustment_3p)
        if (total_abs_length_adj > 0):
            isoform_length_adjustments.append( (total_abs_length_adj, isoform.getTranscriptIDInSrcDB(), first_exon_length_adjustment_5p, last_exon_length_adjustment_3p) )

    isoform_length_adjustments.sort(key=itemgetter(0), reverse=True)
    print("Isoform length adjustments:", file=sys.stderr)
    for tup in isoform_length_adjustments:
        print("%d\t%s\t%d\t%d" % tup, file=sys.stderr)


def groupIsoformsInPhases(isoform_models, annotator, tempdir):
    #isoform_models, all_intersect_results, isoform_instances = cPickle.load(open("pre_grouping.pkl",'rb'))
    isoform_instances, all_intersect_results = intersectIsoforms(isoform_models, tempdir, use_pickled=False)
    #cPickle.dump((isoform_models, all_intersect_results, isoform_instances), open("pre_grouping.pkl",'wb'))
    
    isoform_to_group_map = {}
    
    print("INFO: number of intersection results = %d\n" % len(all_intersect_results), file=sys.stderr)

    phase_num = 1
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    formSingletonGroups(all_intersect_results, isoform_instances, isoform_to_group_map, phase_num)

    phase_num = 2
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    required_src_dbs = set(["GENCODE", "RefSeq"])
    # Group overlapping pairs of isoforms that have same HGNC ID
    remaining_intersect_results = groupIsoformsByOfficialLocus(all_intersect_results, isoform_instances, required_src_dbs, isoform_to_group_map, phase_num)
    print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results), file=sys.stderr)

    phase_num = 3
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    improved = True
    while (improved):
        start_num = len(remaining_intersect_results)
        # Add an isoform to a group only if it is identical to an isoform in that group
        remaining_intersect_results = groupIsoformsWithIdenticalInGroup(remaining_intersect_results, isoform_instances, isoform_to_group_map)
        print("\tINFO: intermediate result %da, remaining number of intersection results = %d" % (phase_num, len(remaining_intersect_results)), file=sys.stderr)
            
        # Add an isoform to a group only if the group contains isoform from same database and locus
        remaining_intersect_results = groupIsoformsSameHGNC(remaining_intersect_results, isoform_instances, isoform_to_group_map)
        print("\tINFO: intermediate result %db, remaining number of intersection results = %d" % (phase_num, len(remaining_intersect_results)), file=sys.stderr)

        improved = len(remaining_intersect_results) < start_num
    print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results), file=sys.stderr)

    #phase_num = 
    #print("INFO: grouping isoforms, Phase %d" % phase_num
    #remaining_intersect_results = groupIsoformsSeededByIdenticalNonSelf(remaining_intersect_results, isoform_instances, isoform_to_group_map, phase_num)
    #print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results)

    phase_num = 4
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    remaining_intersect_results = groupIsoformsNotOverlappingAnyGroup3(remaining_intersect_results, isoform_instances, isoform_to_group_map, phase_num)
    print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results), file=sys.stderr)

    phase_num = 5
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    improved = True
    while (improved):
        start_num = len(remaining_intersect_results)

        # Add an isoform to a group only if it is identical to an isoform in that group
        remaining_intersect_results = groupIsoformsWithIdenticalInGroup(remaining_intersect_results, isoform_instances, isoform_to_group_map)
        print("\tINFO: intermediate result %da, remaining number of intersection results = %d" % (phase_num, len(remaining_intersect_results)), file=sys.stderr)

        # Add an isoform to a group only if the group contains isoform from same database and locus
        remaining_intersect_results = groupIsoformsSameHGNC(remaining_intersect_results, isoform_instances, isoform_to_group_map)
        print("\tINFO: intermediate result %db, remaining number of intersection results = %d" % (phase_num, len(remaining_intersect_results)), file=sys.stderr)
            
        improved = len(remaining_intersect_results) < start_num
    print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results), file=sys.stderr)

    #phase_num = 6
    #print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    #remaining_intersect_results = groupSingleExonIsoformsWithOverlappingInGroup(remaining_intersect_results, isoform_instances, isoform_to_group_map)
    #print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results)

    phase_num = 6
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    remaining_intersect_results = formFusionTranscriptGroups(remaining_intersect_results, isoform_instances, isoform_to_group_map, phase_num)
    print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results), file=sys.stderr)
    
    phase_num = 7
    print("INFO: grouping isoforms, Phase %d" % phase_num, file=sys.stderr)
    remaining_intersect_results = groupIsoformsWithMostSimilar2(remaining_intersect_results, isoform_instances, isoform_to_group_map, phase_num)
    print("INFO: remaining number of intersection results = %d\n" % len(remaining_intersect_results), file=sys.stderr)

    assert (len(remaining_intersect_results)==0)
    print("INFO: Finished grouping isoforms", file=sys.stderr)

    return set(isoform_to_group_map.values()) # (isoform_models, ) TEMPORARY RETURN isoform_models


def intersectIsoforms(isoform_models, tempdir, use_pickled=False):
    intersect_results = []
    isoform_instances = {}

    already_seen = set()
    num_already_seen = 0
    
    if (use_pickled):
        isoform_instances, intersect_results = cPickle.load(open("chr19_intersect_results.pkl", "rb"))
    else:
        isoforms_bed = "%s/isoforms_phase1.bed"  % tempdir
        op_bed = open(isoforms_bed, 'w')
        for isoform in isoform_models:
            bed_fields = isoform.getAsBED12Fields("instance_address")
            isoform_instances[bed_fields[3]] = isoform
            op_bed.write("%s\n" % "\t".join(bed_fields))
        op_bed.close()

        print("INFO: intersecting all isoforms", file=sys.stderr)
        pbt_isoforms = pybedtools.BedTool(isoforms_bed)
        for line in pbt_isoforms.intersect(pbt_isoforms, split=True, s=True, wo=True, stream=True):
            isoform1_frac_in_common_sj, isoform2_frac_in_common_sj, sum_frac_of_all_sj_pos = 0, 0, 0
            bp_overlap = float(line.fields[24])
            isoform1 = isoform_instances[line.fields[3]]
            isoform2 = isoform_instances[line.fields[15]]

            if (isoform1 != isoform2 and (isoform1,isoform2) not in already_seen):
                already_seen.add( (isoform2,isoform1) )
                # Exon positions in common
                isoform1_frac_shared_exon = bp_overlap/float(isoform1.getLengthAllExons())
                isoform2_frac_shared_exon = bp_overlap/float(isoform2.getLengthAllExons())
                sum_frac_overlaps = isoform1_frac_shared_exon + isoform2_frac_shared_exon
                num_shared_exon_termini = len(isoform1.getExonTerminiPositions() & isoform2.getExonTerminiPositions())

                # Splice junction positions in common
                isoform1_sj_pos = isoform1.getSpliceJunctionPositions()
                isoform2_sj_pos = isoform2.getSpliceJunctionPositions()

                if (len(isoform1_sj_pos & isoform2_sj_pos) > 0):
                    all_sj_pos = isoform1_sj_pos.union(isoform2_sj_pos)
                    isoform1_frac_of_all_sj = 0.0 if (len(isoform1_sj_pos)==0) else float(len(all_sj_pos & isoform1_sj_pos))/float(len(all_sj_pos))
                    isoform2_frac_of_all_sj = 0.0 if (len(isoform2_sj_pos)==0) else float(len(all_sj_pos & isoform2_sj_pos))/float(len(all_sj_pos))
                    sum_frac_of_all_sj_pos = isoform1_frac_of_all_sj + isoform2_frac_of_all_sj

                    isoform1_frac_in_common_sj = 0.0 if (len(isoform1_sj_pos)==0) else float(len(all_sj_pos & isoform1_sj_pos))/float(len(isoform1_sj_pos))
                    isoform2_frac_in_common_sj = 0.0 if (len(isoform2_sj_pos)==0) else float(len(all_sj_pos & isoform2_sj_pos))/float(len(isoform2_sj_pos))
                intersect_results.append( (line.fields[3], line.fields[15], sum_frac_overlaps, num_shared_exon_termini, sum_frac_of_all_sj_pos,
                                           isoform1_frac_shared_exon, isoform2_frac_shared_exon, isoform1_frac_in_common_sj, isoform2_frac_in_common_sj) )

            elif ((isoform1,isoform2) in already_seen):
                num_already_seen += 1 

        intersect_results.sort(key=itemgetter(4,3,2), reverse=True)

        # Cleanup
        os.remove(isoforms_bed)
        
    return (isoform_instances, intersect_results)


def groupIdenticalIsoforms(intersect_results, isoform_instances, isoform_to_group_map):
    unused_intersect_results = []
    G = nx.Graph()

    equiv_isoforms = defaultdict(set)
    for tup in intersect_results:
        if (tup[2] == 2.0 or tup[4] == 2.0): # same single-exon isoforms, or multi-exon isoforms with same splice junctions
            isoform1 = isoform_instances[tup[0]]
            isoform2 = isoform_instances[tup[1]]
            G.add_edge(isoform1, isoform2)
            #equiv_isoforms[isoform1].add(isoform2)
            #equiv_isoforms[isoform2].add(isoform1)
        else:
            unused_intersect_results.append(tup)

    for clq in nx.find_cliques(G):
        ig = IsoformGroup()
        for isoform in clq:
            assert (not isoform_to_group_map.has_key(isoform))
            ig.addIsoform(isoform)
            isoform_to_group_map[isoform] = ig

    return unused_intersect_results


def formSingletonGroups(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    '''mRNAs that do not overlap any other mRNAs (on the same strand) form their own group.
    Note that these isoforms are not in the intersection results, so must be identified by their absence.'''
    non_isolated_isoforms = set()
    for tup in intersect_results:
        non_isolated_isoforms.add(tup[0])
        non_isolated_isoforms.add(tup[1])

    for signature, isoform in isoform_instances.items():
        if (signature not in non_isolated_isoforms):
            ig = IsoformGroup(phase_num)
            ig.addIsoform(isoform)
            assert (not isoform_to_group_map.has_key(isoform))
            isoform_to_group_map[isoform] = ig


def groupIsoformsByOfficialLocus(intersect_results, isoform_instances, required_src_dbs, isoform_to_group_map, phase_num):
    unused_intersect_results = []
    G = nx.Graph() # Isoform in different genomic locations and strands can have same HGNC ID, so need to do grouping with graph
    
    for tup in intersect_results:
        tup_used = True
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]

        src_db1, locus1 = isoform1.getDatabaseAndLocus()
        src_db2, locus2 = isoform2.getDatabaseAndLocus()

        isoform1_hgnc_id = isoform1.getHGNCID()
        isoform2_hgnc_id = isoform2.getHGNCID()

        num_exons_1 = isoform1.getNumberOfExons()
        num_exons_2 = isoform2.getNumberOfExons()
        both_single_exon = num_exons_1 == num_exons_2 == 1

        if (src_db1 in required_src_dbs):
            G.add_node(isoform1)

        if (src_db2 in required_src_dbs):
            G.add_node(isoform2)

        if (isoform1_hgnc_id == isoform2_hgnc_id != None):
            G.add_edge(isoform1, isoform2)

        elif (src_db1 == src_db2 and src_db1 in required_src_dbs and locus1 == locus2 != None):
            G.add_edge(isoform1, isoform2)

        elif (src_db1 in required_src_dbs and src_db2 in required_src_dbs and
              (tup[2]==2.0 or tup[4]==2.0 or (both_single_exon and (isoform1_hgnc_id == None or isoform2_hgnc_id == None)))): # Captures single-exon transcript overlaps that aren't quite identical
            G.add_edge(isoform1, isoform2)

        #elif (src_db1 in required_src_dbs and src_db2 in required_src_dbs and both_single_exon):
        #    #if (isoform1_hgnc_id != isoform2_hgnc_id and isoform1_hgnc_id != None and isoform2_hgnc_id != None):
        #        #pdb.set_trace()
        #    G.add_edge(isoform1, isoform2)
                
        else:
            tup_used = False
            
        if (not tup_used):
            unused_intersect_results.append(tup)

    for cc in nx.connected_components(G):
        ig = IsoformGroup(phase_num)
        for isoform in cc:
            try:
                assert (not isoform_to_group_map.has_key(isoform))
            except AssertionError as ae:
                pdb.set_trace()
            found_HGNC_conflict = ig.addIsoform(isoform)
            #if (found_HGNC_conflict):
            #    pdb.set_trace()
            isoform_to_group_map[isoform] = ig

    return unused_intersect_results

    
def groupIsoformsWithIdenticalInGroup(intersect_results, isoform_instances, isoform_to_group_map):
    unused_intersect_results = []
    allow_hgnc_id_conflicts = True

    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]

        num_exons_1 = isoform1.getNumberOfExons()
        num_exons_2 = isoform2.getNumberOfExons()
        both_single_exon = num_exons_1 == num_exons_2 == 1

        # tup[7]==1.0 or tup[8]==1.0 or 
        if (tup[2] == 2.0 or tup[4] == 2.0 or (both_single_exon and tup[2]>1.95)): # Captures single-exon transcript overlaps that aren't quite identical
            isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
            isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

            if ( xor(isoform1_group == None, isoform2_group == None) ):
                if (isoform1_group == None):
                    found_HGNC_conflict = isoform2_group.addIsoform(isoform1, allow_hgnc_id_conflicts)
                    if (found_HGNC_conflict):
                        pdb.set_trace()
                    isoform_to_group_map[isoform1] = isoform2_group
                else:
                    found_HGNC_conflict = isoform1_group.addIsoform(isoform2, allow_hgnc_id_conflicts)
                    if (found_HGNC_conflict):
                        pdb.set_trace()
                    isoform_to_group_map[isoform2] = isoform1_group

            elif (isoform1_group != None and isoform2_group != None and isoform1_group != isoform2_group):
                group1_symbol = isoform1_group.getHGNCSymbol()
                group2_symbol = isoform2_group.getHGNCSymbol()
                print("\tINFO: recording conflict due to isoforms overlap between isoform groups %s and %s" % (group1_symbol, group2_symbol), file=sys.stderr)
                isoform1_group.recordHGNCConflictWith( isoform2_group )
                isoform2_group.recordHGNCConflictWith( isoform1_group )
            else:
                unused_intersect_results.append(tup)
        else:
            unused_intersect_results.append(tup)

    #unused_intersect_results = filter(lambda tup: tup[0] != tup[1] or isoform_to_group_map[isoform_instances[tup[0]]] == None, unused_intersect_results)
    return unused_intersect_results


def groupIsoformsSameHGNC(intersect_results, isoform_instances, isoform_to_group_map):
    unused_intersect_results = []

    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None

        isoform2 = isoform_instances[tup[1]]
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if (isoform1.getHGNCID() == isoform2.getHGNCID() != None):
            if (isoform1_group == isoform2_group == None):
                unused_intersect_results.append(tup)
            elif (isoform1_group == None):
                found_HGNC_conflict = isoform2_group.addIsoform(isoform1)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                isoform_to_group_map[isoform1] = isoform2_group
            elif (isoform2_group == None):
                found_HGNC_conflict = isoform1_group.addIsoform(isoform2)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                isoform_to_group_map[isoform2] = isoform1_group
            else:
                try:
                    assert (isoform1_group==isoform2_group), "Isoform groups not the same"
                except AssertionError:
                    pdb.set_trace()
        else:
            unused_intersect_results.append(tup)

    return unused_intersect_results


def groupIsoformsSameDBAndLocus(intersect_results, isoform_instances, isoform_to_group_map):
    '''SameDB is used here to mean same locus of same database or same standardized locus symbol or HGNC ID.'''
    unused_intersect_results = []

    # Determine the non-None isoform groups that each isoform matches to (by having a common DB+locus).
    # If an isoform matches to more than one group, then it is a fusion transcript between multiple canonical
    # loci and is destined to become its own group.
    group_matches_per_isoform = defaultdict(set)
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None

        isoform2 = isoform_instances[tup[1]]
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if (isoform2_group != None and (isoform2_group.containsIsoformFromSameDBAndLocusAs(isoform1) or
                                        isoform1.getLocusSymbol() == isoform2.getLocusSymbol() != None or
                                        isoform1.getHGNCID() == isoform2.getHGNCID() != None)):
            group_matches_per_isoform[isoform1].add(isoform2_group)

        if (isoform1_group != None and (isoform1_group.containsIsoformFromSameDBAndLocusAs(isoform2) or
                                        isoform1.getLocusSymbol() == isoform2.getLocusSymbol() != None or
                                        isoform1.getHGNCID() == isoform2.getHGNCID() != None)):
            group_matches_per_isoform[isoform2].add(isoform1_group)

    for tup in intersect_results:
        tup_used = False
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]

        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if ( xor(isoform1_group == None, isoform2_group == None) ):
            if (isoform2_group != None and (isoform1.getHGNCID() == isoform2.getHGNCID() != None or
                                            (len(group_matches_per_isoform[isoform1])==1 and isoform2_group in group_matches_per_isoform[isoform1]))):
                tup_used = True
                found_HGNC_conflict = isoform2_group.addIsoform(isoform1)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                isoform_to_group_map[isoform1] = isoform2_group

            elif (isoform1_group != None and (isoform1.getHGNCID() == isoform2.getHGNCID() != None or
                                              (len(group_matches_per_isoform[isoform2])==1 and isoform1_group in group_matches_per_isoform[isoform2]))):
                tup_used = True
                found_HGNC_conflict = isoform1_group.addIsoform(isoform2)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                isoform_to_group_map[isoform2] = isoform1_group

        elif (isoform1_group == isoform2_group != None):
            #print >> sys.stderr, isoform1_group.hgnc_id
            #print >> sys.stderr, isoform1_group.db_genes
            #print >> sys.stderr, isoform2_group.hgnc_id
            #print >> sys.stderr, isoform2_group.db_genes
            #pdb.set_trace()
            #print >> sys.stderr, "-----------------------------"
            tup_used = True

        if (not tup_used):
            unused_intersect_results.append(tup)

    #unused_intersect_results = filter(lambda tup: tup[0] != tup[1] or isoform_to_group_map[isoform_instances[tup[0]]] == None, unused_intersect_results)
    return unused_intersect_results


def groupSingleExonIsoformsWithOverlappingInGroup(intersect_results, isoform_instances, isoform_to_group_map):
    unused_intersect_results = []
    allow_hgnc_id_conflicts = True

    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]

        num_exons_1 = isoform1.getNumberOfExons()
        num_exons_2 = isoform2.getNumberOfExons()

        if (num_exons_1 == num_exons_2 == 1):
            isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
            isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

            if ( xor(isoform1_group == None, isoform2_group == None) ):
                if (isoform1_group == None):
                    found_HGNC_conflict = isoform2_group.addIsoform(isoform1, allow_hgnc_id_conflicts)
                    if (found_HGNC_conflict):
                        pdb.set_trace()
                    isoform_to_group_map[isoform1] = isoform2_group
                else:
                    found_HGNC_conflict = isoform1_group.addIsoform(isoform2, allow_hgnc_id_conflicts)
                    if (found_HGNC_conflict):
                        pdb.set_trace()
                    isoform_to_group_map[isoform2] = isoform1_group

            #elif (isoform1_group != None and isoform2_group != None and isoform1_group != isoform2_group):
            #    isoform1_hgnc_id = isoform1.getHGNCID()
            #    isoform2_hgnc_id = isoform2.getHGNCID()
            #    if (isoform1_hgnc_id != None and isoform2_hgnc_id != None):
            #        pdb.set_trace()
            #        isoform1_group.recordHGNCConflictWith( isoform2_group )
            #        isoform2_group.recordHGNCConflictWith( isoform1_group )
            else:
                unused_intersect_results.append(tup)
        else:
            unused_intersect_results.append(tup)

    #unused_intersect_results = filter(lambda tup: tup[0] != tup[1] or isoform_to_group_map[isoform_instances[tup[0]]] == None, unused_intersect_results)
    return unused_intersect_results


def formFusionTranscriptGroups(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []

    # For each isoform not in a group, determine the non-None isoform groups that is overlaps
    group_matches_per_isoform = defaultdict(set)
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
        isoform2 = isoform_instances[tup[1]]
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if (isoform1_group == None and isoform2_group != None and tup[3] > 0):
            group_matches_per_isoform[isoform1].add(isoform2_group)

        if (isoform1_group != None and isoform2_group == None and tup[3] > 0):
            group_matches_per_isoform[isoform2].add(isoform1_group)

            
    # By definition, fusion transcripts overlap multiple isoform groups.
    # Those transcripts that only overlap one shouldn't be handled in this method.
    isoforms_to_ignore = set()
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
        isoform2 = isoform_instances[tup[1]]
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if (isoform1_group == None and isoform2_group != None and tup[3] > 0):
            if (len(group_matches_per_isoform[isoform1]) == 1):
                isoforms_to_ignore.add(isoform1)
                
        if (isoform2_group == None and isoform1_group != None and tup[3] > 0):
            if (len(group_matches_per_isoform[isoform2]) == 1):
                isoforms_to_ignore.add(isoform2)
            
    for isoform in isoforms_to_ignore:
        del group_matches_per_isoform[isoform]            


    isoforms_w_same_group_match = defaultdict(list)
    for isoform, isoform_group_set in group_matches_per_isoform.items():
        isoforms_w_same_group_match[frozenset(isoform_group_set)].append(isoform)
        #isoform.setAsFusionTranscript()

    G = nx.Graph()
    G.add_nodes_from(isoforms_w_same_group_match.keys())
    for ig_set_A, ig_set_B in combinations(isoforms_w_same_group_match.keys(),2):
        if (ig_set_A.issubset(ig_set_B) or ig_set_B.issubset(ig_set_A)): # So that isoforms overlapping groups A,B,C are merged with those overlapping groups A,B
            G.add_edge(ig_set_A, ig_set_B)

    for cc in nx.connected_components(G):
        # Only form new isoform group if cc contains RefSeq/HGNC/GENCODE isoforms.
        new_ig_okay = False
        for isoform in chain.from_iterable(map(lambda k: isoforms_w_same_group_match[k], cc)):
            refseq_locus = isoform.getLocusSymbolFromDatabase("RefSeq")
            gencode_locus = isoform.getLocusSymbolFromDatabase("GENCODE")
            hgnc_id = isoform.getHGNCID()
            if (refseq_locus != None or gencode_locus != None or hgnc_id != None):
                new_ig_okay = True
                break
            
        if (new_ig_okay):
            ig = IsoformGroup(phase_num)
            for isoform in chain.from_iterable(map(lambda k: isoforms_w_same_group_match[k], cc)):
                found_HGNC_conflict = ig.addIsoform(isoform)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                assert (not isoform_to_group_map.has_key(isoform)), "Isoform %s already in a group" % isoform.getTranscriptIDInSrcDB()
                isoform_to_group_map[isoform] = ig

    # Keep the intersect results that still have isoforms that are not associated with an isoform group
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
        isoform2 = isoform_instances[tup[1]]
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if (isoform1_group == None or isoform2_group == None):
            unused_intersect_results.append(tup)
    unused_intersect_results.sort(key=itemgetter(3,4,2), reverse=True)

    return unused_intersect_results


def groupIsoformsSeededByIdenticalNonSelf(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []

    # Collect isoforms with identical sequence and/or splice junctions
    equiv_isoforms = defaultdict(set)
    for tup in intersect_results:
        try:
            assert (tup[0] != tup[1]) # TODO Should be able to remove
        except AssertionError:
            pdb.set_trace()
            
        if (tup[2] == 2.0 or tup[4] == 2.0): #  and tup[0] != tup[1]): REMOVE
            isoform1 = isoform_instances[tup[0]]
            isoform2 = isoform_instances[tup[1]]
            equiv_isoforms[isoform1].add(isoform2)
            equiv_isoforms[isoform2].add(isoform1)
            
    for tup in intersect_results:
        if (tup[2] == 2.0 or tup[4] == 2.0): #  and tup[0] != tup[1] TODO REMOVE
            isoform1 = isoform_instances[tup[0]]
            isoform2 = isoform_instances[tup[1]]
            
            isoform1_group = isoform_to_group_map[isoform1]
            isoform2_group = isoform_to_group_map[isoform2]

            add_isoform1, add_isoform2 = False, False
            ig = None

            if (isoform1_group != None and isoform2_group == None):
                add_isoform2 = True
                ig = isoform1_group
            elif (isoform1_group == None and isoform2_group != None):
                add_isoform1 = True
                ig = isoform2_group
            elif (isoform1_group == None and isoform2_group == None):
                ig = IsoformGroup(phase_num)
                add_isoform1 = True
                # isoform2 gets added below because it is equivalent to isoform1.
            else:
                assert (isoform1_group == isoform2_group != None)
                #isoform1_group.addIntersectResult(tup)
                
            if (add_isoform1):
                assert (not isoform_to_group_map.has_key(isoform1))
                found_HGNC_conflict = ig.addIsoform(isoform1)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()

                isoform_to_group_map[isoform1] = ig
                for equiv_isoform in equiv_isoforms[isoform1]:
                    if (not ig.hasIsoform(equiv_isoform)):
                        assert (not isoform_to_group_map.has_key(equiv_isoform))
                        found_HGNC_conflict = ig.addIsoform(equiv_isoform)
                        #if (found_HGNC_conflict):
                        #    pdb.set_trace()

                        assert (not isoform_to_group_map.has_key(equiv_isoform))
                        isoform_to_group_map[equiv_isoform] = ig
                    else:
                        assert(isoform_to_group_map[equiv_isoform] == ig)

            if (add_isoform2):
                assert (not isoform_to_group_map.has_key(isoform2))
                found_HGNC_conflict = ig.addIsoform(isoform2)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()

                isoform_to_group_map[isoform2] = ig
                for equiv_isoform in equiv_isoforms[isoform2]:
                    if (not ig.hasIsoform(equiv_isoform)):
                        assert (not isoform_to_group_map.has_key(equiv_isoform))
                        found_HGNC_conflict = ig.addIsoform(equiv_isoform)
                        #if (found_HGNC_conflict):
                        #    pdb.set_trace()

                        isoform_to_group_map[equiv_isoform] = ig
                    else:
                        assert(isoform_to_group_map[equiv_isoform] == ig)

        elif (not isoform_to_group_map.has_key(isoform_instances[tup[0]])):
            unused_intersect_results.append(tup)

    return unused_intersect_results


def groupIsoformsNotOverlappingAnyGroup(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []
    
    intersect_results_grouped_by_isoform1 = OrderedDict()
    for tup in intersect_results:
        isoform1_name = tup[0]
        if (not isoform_to_group_map.has_key(isoform_instances[isoform1_name])):
            if (not intersect_results_grouped_by_isoform1.has_key(isoform1_name)):
                intersect_results_grouped_by_isoform1[isoform1_name] = []
            intersect_results_grouped_by_isoform1[isoform1_name].append(tup)

    for isoform1_name, intersect_results_for_isoform1 in intersect_results_grouped_by_isoform1.items():
        overlapped_group_found = False
        for tup in intersect_results_for_isoform1:
            isoform2 = isoform_instances[tup[1]]
            if (isoform_to_group_map.has_key(isoform2)):
                overlapped_group_found = True
                break

        if (not overlapped_group_found):
            isoform1 = isoform_instances[isoform1_name]
            new_ig = IsoformGroup(phase_num)
            assert (not isoform_to_group_map.has_key(isoform1))
            found_HGNC_conflict = new_ig.addIsoform(isoform1)
            if (found_HGNC_conflict):
                pdb.set_trace()

            isoform_to_group_map[isoform1] = new_ig
            for tup in intersect_results_for_isoform1:
                isoform2 = isoform_instances[tup[1]]

                num_exons_1 = isoform1.getNumberOfExons()
                num_exons_2 = isoform2.getNumberOfExons()
                both_single_exon = num_exons_1 == num_exons_2 == 1

                if (tup[2] == 2.0 or tup[4] == 2.0 or isoform2.isFromSameDBAndLocusAs(isoform1) or (both_single_exon and tup[2]>1.95)):
                    assert (not isoform_to_group_map.has_key(isoform2))
                    found_HGNC_conflict = new_ig.addIsoform(isoform2)
                    if (found_HGNC_conflict):
                        pdb.set_trace()
                    isoform_to_group_map[isoform2] = new_ig
                else:
                    unused_intersect_results.append(tup)
        else:
            unused_intersect_results.extend(intersect_results_for_isoform1)

    #unused_intersect_results = filter(lambda tup: tup[0] != tup[1] or isoform_to_group_map[isoform_instances[tup[0]]] == None, unused_intersect_results)    
    unused_intersect_results = sorted(unused_intersect_results, key=itemgetter(4,2,3), reverse=True)
    return unused_intersect_results


def groupIsoformsNotOverlappingAnyGroup2(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []

    isoforms_that_overlap_a_group = set()
    
    # Make graph that connects compatibly-overlapping isoforms
    G = nx.Graph()
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]

        num_exons_1 = isoform1.getNumberOfExons()
        num_exons_2 = isoform2.getNumberOfExons()
        both_single_exon = num_exons_1 == num_exons_2 == 1
        both_multi_exon = num_exons_1 == num_exons_2 > 1

        if (not isoform_to_group_map.has_key(isoform1) and not isoform_to_group_map.has_key(isoform2)):
            if (tup[3]>0 or both_single_exon or isoform1.getLocusSymbol() == isoform2.getLocusSymbol()):
                G.add_edge(isoform1,isoform2)
            else:
                G.add_nodes_from([isoform1,isoform2])
        elif (not isoform_to_group_map.has_key(isoform1) and isoform_to_group_map.has_key(isoform2)):
            if (tup[3]==0):
                G.add_node(isoform1)
            else:
                isoforms_that_overlap_a_group.add(isoform1)
        elif (isoform_to_group_map.has_key(isoform1) and not isoform_to_group_map.has_key(isoform2)):
            if (tup[3]==0):
                G.add_node(isoform2)
            else:
                isoforms_that_overlap_a_group.add(isoform2)

    # Get connected components that represent new groups of isoforms that don't overlap existing groups.
    G.remove_nodes_from(isoforms_that_overlap_a_group)
    for cc in nx.connected_components(G):
        new_ig = IsoformGroup(phase_num)
        for isoform in cc:
            new_ig.addIsoform(isoform)
            assert (not isoform_to_group_map.has_key(isoform))
            isoform_to_group_map[isoform] = new_ig

    # Intersection results that have a corresponding edge in G have been used. Get the remaining unused intersection results
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]
        if (not G.has_edge(isoform1, isoform2)):
            unused_intersect_results.append(tup)
    unused_intersect_results.sort(key=itemgetter(4,2,3), reverse=True)

    return unused_intersect_results


def groupIsoformsNotOverlappingAnyGroup3(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []

    # Make graph that connects compatibly-overlapping isoforms
    G = nx.Graph()
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]

        #num_exons_1 = isoform1.getNumberOfExons()
        #num_exons_2 = isoform2.getNumberOfExons()
        #both_single_exon = num_exons_1 == num_exons_2 == 1
        #both_multi_exon = num_exons_1 == num_exons_2 > 1

        #if (tup[3]>0 or not both_multi_exon or isoform1.getLocusSymbol() == isoform2.getLocusSymbol()):
        G.add_edge(isoform1,isoform2)            

    # Get connected components that represent new groups of isoforms that don't overlap existing groups.
    for cc in nx.connected_components(G):
        if (all(map(lambda isoform: isoform not in isoform_to_group_map, cc))): # If no isoform in the connected component is in a group...
            new_ig = IsoformGroup(phase_num)
            for isoform in cc:
                new_ig.addIsoform(isoform)
                assert (not isoform_to_group_map.has_key(isoform))
                isoform_to_group_map[isoform] = new_ig

            for isoform1, isoform2 in combinations(cc,2):
                if (G.has_edge(isoform1,isoform2)):
                    G.remove_edge(isoform1,isoform2)
                
    # Intersection results that have a corresponding edge in G have been used. Get the remaining unused intersection results
    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]
        if (G.has_edge(isoform1, isoform2)):
            unused_intersect_results.append(tup)
    unused_intersect_results.sort(key=itemgetter(4,2,3), reverse=True)

    return unused_intersect_results


def groupIsoformsWithMostSimilar(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []

    for tup in intersect_results:
        isoform1 = isoform_instances[tup[0]]
        isoform2 = isoform_instances[tup[1]]
            
        isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
        isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

        if ( xor(isoform1_group == None, isoform2_group == None) ):
            if (isoform1_group == None):
                found_HGNC_conflict = isoform2_group.addIsoform(isoform1)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                isoform_to_group_map[isoform1] = isoform2_group
            else:
                found_HGNC_conflict = isoform1_group.addIsoform(isoform2)
                #if (found_HGNC_conflict):
                #    pdb.set_trace()
                isoform_to_group_map[isoform2] = isoform1_group
        else:
            unused_intersect_results.append(tup)

        #elif (isoform1_group == None and isoform2_group == None):
        #    new_ig = IsoformGroup(phase_num)
        #    new_ig.addIsoform(isoform1)
        #    isoform_to_group_map[isoform1] = new_ig
        #    new_ig.addIsoform(isoform2)
        #    isoform_to_group_map[isoform2] = new_ig
        #
        #    unused_intersect_results.append(tup)

        #elif (isoform1_group != isoform2_group):
        #    pdb.set_trace()
        #    print isoform1.ID, isoform2.ID
            
    # All non-identical overlapping isoforms that were both not part of a group should be resolved at this point.
    unused_intersect_results = filter(lambda tup: not (isoform_to_group_map.has_key(isoform_instances[tup[0]]) and
                                      isoform_to_group_map.has_key(isoform_instances[tup[1]])), unused_intersect_results)    
    return unused_intersect_results


def groupIsoformsWithMostSimilar2(intersect_results, isoform_instances, isoform_to_group_map, phase_num):
    unused_intersect_results = []

    num_begin = len(intersect_results)
    num_end = -1
    
    while (num_begin > num_end):
        unused_intersect_results = []
        num_begin = len(intersect_results)
        for tup in intersect_results:
            isoform1 = isoform_instances[tup[0]]
            isoform2 = isoform_instances[tup[1]]
            
            isoform1_group = isoform_to_group_map[isoform1] if isoform_to_group_map.has_key(isoform1) else None
            isoform2_group = isoform_to_group_map[isoform2] if isoform_to_group_map.has_key(isoform2) else None

            if ( xor(isoform1_group == None, isoform2_group == None) ):
                if (isoform1_group == None):
                    found_HGNC_conflict = isoform2_group.addIsoform(isoform1)
                    #if (found_HGNC_conflict):
                    #    pdb.set_trace()
                    isoform_to_group_map[isoform1] = isoform2_group
                else:
                    found_HGNC_conflict = isoform1_group.addIsoform(isoform2)
                    #if (found_HGNC_conflict):
                    #    pdb.set_trace()
                    isoform_to_group_map[isoform2] = isoform1_group
            else:
                unused_intersect_results.append(tup)

            #elif (isoform1_group == None and isoform2_group == None):
            #    new_ig = IsoformGroup(phase_num)
            #    new_ig.addIsoform(isoform1)
            #    isoform_to_group_map[isoform1] = new_ig
            #    new_ig.addIsoform(isoform2)
            #    isoform_to_group_map[isoform2] = new_ig
            #
            #    unused_intersect_results.append(tup)
            #
            #elif (isoform1_group != isoform2_group):
            #    pdb.set_trace()
            #    print isoform1.ID, isoform2.ID
            
        # All non-identical overlapping isoforms that were both not part of a group should be resolved at this point.
        unused_intersect_results = filter(lambda tup: not (isoform_to_group_map.has_key(isoform_instances[tup[0]]) and
                                                           isoform_to_group_map.has_key(isoform_instances[tup[1]])), unused_intersect_results)    
        num_end = len(unused_intersect_results)
        intersect_results = unused_intersect_results[:]
        
    return unused_intersect_results



def debugStandardizePromoterRegions(genome_seq):
    (isoform_models, G, closest_upstream_tss) = cPickle.load(open("chr8.pkl", "rb"))
    print("INFO: standardizing promoter regions", file=sys.stderr)

    for component in nx.connected_components(G):
        if (len(component)>1 and any(map(lambda c: isinstance(c, Exon), component))):
            standardizePromoterRegion(component, closest_upstream_tss, genome_seq, isoform_models)
    print("done", file=sys.stderr)
    return isoform_models


def debugStandardizePolyASites():
    (isoform_models, G) = cPickle.load(open("chr8_pAs.pkl", "rb"))
    print("INFO: standardizing polyA sites", file=sys.stderr)

    for pAs_and_exons in nx.connected_components(G):
        if (len(pAs_and_exons)>1 and any(map(lambda c: isinstance(c, Exon) and c.isThreePrimeComplete(), pAs_and_exons))):
            standardizePolyASites(pAs_and_exons)

    print("done", file=sys.stderr)
    return isoform_models


def consolidateGroupedIsoforms(isoform_groups, genome_seq): # , isoform_models
    for counter, ig in enumerate(isoform_groups):
        #print >> sys.stderr, "\rig %d" % counter,
        ig.consolidateIsoforms(genome_seq) # isoform_models
    #print >> sys.stderr, "\n"


def dumpAbsorptionInfo(all_isoform_models, ordered_chromosomes, filename):
    op = open(filename, 'w')
    for chrom in ordered_chromosomes:
        for isoform in all_isoform_models[chrom]:

            ID = isoform.getTranscriptIDInSrcDB()
            ID = str(ID) # In case 'None'
            CGID = isoform.getCGDBName() if isoform.hasCGDBName() else None
            absorbed = isoform.getAbsorbed()
            absorbed_by = isoform.getAbsorbedBy()

            if (len(absorbed_by) == 0):
                assert (CGID != None)

            if (len(absorbed) == 0 and len(absorbed_by) == 0):
                op.write("%s <- %s\n" % (CGID, ID))

            # Write absorption into this isoform
            for other_isoform in absorbed:
                other_ID = other_isoform.getTranscriptIDInSrcDB()        
                assert (not other_isoform.hasCGDBName())
                op.write("%s (%s) <- %s\n" % (CGID, ID, other_ID))

            # Write which isoform(s) this isoform was absorbed into (should be redundant info to absorbed)
            #for other_isoform in absorbed_by:
            #    other_absorbed_by = other_isoform.getAbsorbedBy()
            #    if (len(other_absorbed_by) == 0):
            #        other_ID = other_isoform.getTranscriptIDInSrcDB()
            #        assert (other_isoform.hasCGDBName())
            #        other_CGID = other_isoform.getCGDBName()
            #        op.write("%s (%s) <= %s\n" % (other_CGID, other_ID, ID))            
    op.close()


def verifyIsoformGroups(isoform_groups, all_isoform_models):
    # Make certain that each isoform is accounted for in some isoform group
    all_ig_isoforms = set()
    for ig in isoform_groups:
        ig_isoforms = ig.getIsoforms()
        all_ig_isoforms.update(ig_isoforms)

    for isoform in all_isoform_models:
        assert (isoform in all_ig_isoforms), "ERROR: %s isoform %s:%s was is not in any IsoformGroup" % (isoform.getDatabaseAndLocus() + (isoform.getTranscriptIDInSrcDB(),))
                              
    # Other checks
    for counter, ig in enumerate(isoform_groups):
        #print >> sys.stderr, "\rig %d" % counter,
        ig.performFinalChecks(all_isoform_models)
    #print >> sys.stderr, "\n"

    
def labelGroupsAndIsoforms(isoform_groups_by_chromosome, ordered_chromosomes):
    group_counter = 1
    
    for chromosome in ordered_chromosomes:
        print("INFO: labelling %s" % chromosome, file=sys.stderr)
        sorted_isoform_groups = sorted(isoform_groups_by_chromosome[chromosome], key=methodcaller("getFivePrimeBoundary"))
        for ig in sorted_isoform_groups:
            ig.setGroupNumber(group_counter)
            ig.setIsoformTSSAndSpliceVariantPolyASiteNumbers()
            for isoform in ig.getIsoforms(not_absorbed=True):
                isoform.setCGDBName(group_counter)
            ig.setLocusSymbol()
            group_counter += 1
            
            
def writeCGDB(isoform_groups_by_chromosome, ordered_chromosomes, chromosome_sizes_tsv, output_base_filename):

    #pkl_file = "%s.pkl" % output_base_filename
    #print >> sys.stderr, "\nINFO: writing %s" % pkl_file
    #cPickle.dump((isoform_groups_by_chromosome, ordered_chromosomes, chromosome_sizes_tsv), open(pkl_file, "wb"))

    chromosome_sizes = {}
    with open(chromosome_sizes_tsv, 'r') as ip:
        header = ip.readline()
        for line in ip:
            chrom_name, chrom_size = line.strip().split("\t")
            chromosome_sizes[chrom_name] = chrom_size
            
    print("INFO: writing output GFF3, BED, and GTF files with basename %s" % output_base_filename, file=sys.stderr)

    op_bed = open("%s.bed" % output_base_filename, 'w')
    op_gff3 = open("%s.gff3" % output_base_filename, 'w')
    op_gtf = open("%s.gtf" % output_base_filename, 'w')

    op_gff3.write("##gff-version 3\n")
    op_gff3.write("##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606\n")

    for chromosome in ordered_chromosomes:
        op_gff3.write("##sequence-region %s 1 %s\n" % (chromosome, chromosome_sizes[chromosome]))

        sorted_isoform_groups = sorted(isoform_groups_by_chromosome[chromosome], key=attrgetter("group_number"))

        for ig in sorted_isoform_groups:
            ig.writeAsBED(op_bed, write_absorbed=False)
            ig.writeAsGFF3(op_gff3, write_absorbed=False)
            ig.writeAsGTF(op_gtf, write_absorbed=False)

    op_bed.close()
    op_gff3.close()
    op_gtf.close()
    
    
if (__name__ == "__main__"):
    tempdir_root, genome_fasta, tss_databases, isoform_model_databases, polyAsite_databases, chromosome_sizes_tsv, \
    hgnc_data_tsv, output_base_filename = sys.argv[1:]
    
    tempdir = "%s/createCGDB_%s_%d" % (tempdir_root, os.getlogin(), os.getpid())
    os.mkdir(tempdir)
    pybedtools.set_tempdir(tempdir)

    genome_seq = Fasta(genome_fasta)

    annotator = RNAIsoformAnnotator(hgnc_data_tsv)
    all_isoform_models = readIsoformModels(isoform_model_databases, annotator)
    all_TSSs = readTSSs(tss_databases)
    all_polyAsites = readPolyASites(polyAsite_databases)
    ordered_chromosomes = orderChromosomes(all_TSSs, all_isoform_models)
            
    isoform_groups_by_chromosome = {}
    ordered_chromosomes = ["chr19"]
    for chromosome in ordered_chromosomes:
        print("\nINFO: processing %s" % chromosome, file=sys.stderr)
        isoform_models = all_isoform_models[chromosome]

        print("INFO: grouping and standardizing TSS", file=sys.stderr)
        TSSs = all_TSSs[chromosome]
        groupTranscriptionStarts(isoform_models, TSSs, genome_seq, tempdir)

        print("INFO: grouping and standardizing polyA sites", file=sys.stderr)
        pAs = all_polyAsites[chromosome]
        groupPolyASites(isoform_models, pAs, tempdir)

        print("INFO: setting final isoform lengths and exon termini", file=sys.stderr)
        for isoform in isoform_models:
            isoform.setLengthAndExonTermini()
            isoform.setUniqueSignature()

        #reportIsoformTerminiLengthAdjustments(isoform_models)

        #isoform_models = cPickle.load(open("isoform_models.pkl",'rb')) # TEMPORARY FOR DEBUGGING
        print("INFO: grouping isoforms", file=sys.stderr)
        isoform_groups = groupIsoformsInPhases(isoform_models, annotator, tempdir)
        #cPickle.dump(isoform_groups, open("%s_isoform_groups.pkl" % chromosome, "wb"))
        #cPickle.dump((ordered_chromosomes, annotator, isoform_models, isoform_groups), open("%s_after_IG_creation.pkl" % chromosome, "wb"))

        print("INFO: merging intra-group isoforms", file=sys.stderr)
        consolidateGroupedIsoforms(isoform_groups, genome_seq) # , isoform_models
        #isoform_groups = cPickle.load(open("%s_isoform_groups.pkl" % chromosome, "rb"))

        print("INFO: verifying isoform groups", file=sys.stderr)
        verifyIsoformGroups(isoform_groups, isoform_models)
        #cPickle.dump((ordered_chromosomes, annotator, isoform_models, isoform_groups), open("%s_after_IG_consolidation.pkl" % chromosome, "wb"))

        isoform_groups_by_chromosome[chromosome] = isoform_groups
        pybedtools.cleanup()

    labelGroupsAndIsoforms(isoform_groups_by_chromosome, ordered_chromosomes)

    writeCGDB(isoform_groups_by_chromosome, ordered_chromosomes, chromosome_sizes_tsv, output_base_filename)

    filename = "%s_absorption_info.txt" % output_base_filename
    dumpAbsorptionInfo(all_isoform_models, set(ordered_chromosomes), filename)

    shutil.rmtree(tempdir, ignore_errors=True)

    sys.exit(0)
    
