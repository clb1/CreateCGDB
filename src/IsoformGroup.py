import sys
from itertools import product
from operator import itemgetter, methodcaller
from collections import defaultdict
from RNAIsoform_python2 import RNAIsoform
import networkx as nx

import pdb

class IsoformGroup(object):
    def __init__(self, creation_phase):
        self.creation_phase = creation_phase

        self.locus_symbol = None
        self.group_number = None
        self.chromosome = None
        self.strand = None
        self.isoforms = set()
        self.isoforms_IDs = set()
        self.db_genes = {} # There should only be one gene symbol (value) per database (key)

        self.start = 1e30
        self.stop = -1
        
        #self.gene_symbols = set()
        self.hgnc_id = None
        self.hgnc_symbol = None

        self.conflicting_hgnc_ids = []
        self.conflicting_hgnc_symbols = []
        self.conflicting_db_genes = defaultdict(set)

        #self.intersect_results = set()   # TODO: Delete. Not needed.


    def getCreationPhase(self):
        return self.creation_phase


    def addIsoform(self, isoform, known_hgnc_conflict=False):
        hgnc_id = isoform.getHGNCID()
        hgnc_symbol = isoform.getHGNCSymbol()
        isoform_ID = isoform.getTranscriptIDInSrcDB()
        found_HGNC_conflict = self.checkHGNC(hgnc_id, hgnc_symbol, known_hgnc_conflict)

        db, gene_symbol = isoform.getDatabaseAndLocus()
        if (self.db_genes.has_key(db)):
            try:
                assert (gene_symbol == self.db_genes[db])
            except AssertionError:
                self.conflicting_db_genes[db].add(gene_symbol)
        else:
            self.db_genes[db] = gene_symbol
            #self.gene_symbols.add(gene_symbol)
                
        try:
            assert (isoform not in self.isoforms), "Already added isoform %s to IsoformGroup" % isoform_ID
        except AssertionError as ae:
            print >> sys.stderr, ae.message
            pdb.set_trace()

        self.isoforms.add(isoform)
        self.isoforms_IDs.add(isoform_ID)

        chromosome = isoform.getChromosome()
        if (self.chromosome == None):
            self.chromosome = chromosome
        else:
            assert (chromosome == self.chromosome)

        strand = isoform.getStrand()
        if (self.strand == None):
            self.strand = strand
        else:
            assert (strand == self.strand)

        start, stop = isoform.getStartStop()
        if (self.start > start):
            self.start = start

        if (self.stop < stop):
            self.stop = stop
        
        return found_HGNC_conflict


    def addIsoforms(self, isoforms, known_hgnc_conflict=False):
        for isoform in isoforms:
            self.addIsoform(isoform, known_hgnc_conflict)


    #def addIntersectResult(self, intersect_tuple):
    #    self.intersect_results.add(intersect_tuple)


    #def addIntersectResults(self, intersect_tuples):
    #    self.intersect_results.update(intersect_tuples)


    #def mergeWithIsoformGroup(self, isoform_group):
    #    hgnc_id = isoform_group.getHGNCID()
    #    self.checkHGNCID(hgnc_id)

    #    hgnc_symbol = isoform_group.getHGNCSymbol()
    #    self.checkHGNCSymbol(hgnc_symbol)
            
    #    gene_symbols = isoform_group.getAltSymbols()
    #    self.gene_symbols.update(gene_symbols)

    #    new_isoforms = isoform_group.getIsoforms()
    #    self.isoforms.update(new_isoforms)

    def mergeWith(self, other_isoform_group):
        for isoform in other_isoform_group.getIsoforms():
            self.addIsoform(isoform)


    def checkHGNC(self, new_hgnc_id, new_hgnc_symbol, known_hgnc_conflict=False):
        found_hgnc_ID_conflict = False
        if (new_hgnc_id != None and self.hgnc_id != None):
            if (new_hgnc_id != self.hgnc_id):
                if (not known_hgnc_conflict):
                    found_hgnc_ID_conflict = True
                self.conflicting_hgnc_ids.append(new_hgnc_id)
        elif (new_hgnc_id != None and self.hgnc_id == None):
            self.hgnc_id = new_hgnc_id

        found_hgnc_symbol_conflict = False
        if (new_hgnc_symbol != None and self.hgnc_symbol != None):
            if (new_hgnc_symbol != self.hgnc_symbol):
                if (not known_hgnc_conflict):
                    found_hgnc_symbol_conflict = True
                self.conflicting_hgnc_symbols.append(new_hgnc_symbol)
        elif (new_hgnc_symbol != None and self.hgnc_symbol == None):
            self.hgnc_symbol = new_hgnc_symbol

        if (found_hgnc_ID_conflict or found_hgnc_symbol_conflict):
            print >> sys.stderr, "\tWARNING: merging isoform with HGNC Symbol %s (%s) into Isoform Group with HGNC Symbol %s (%s)" % (new_hgnc_symbol, new_hgnc_id, self.hgnc_symbol, self.hgnc_id)
            
        return found_hgnc_ID_conflict or found_hgnc_symbol_conflict


    def recordHGNCConflictWith(self, conflicting_isoform_group):
        conflicting_hgnc_id, conflicting_hgnc_symbol = conflicting_isoform_group.getHGNCIDAndSymbol()

        if (not conflicting_hgnc_id in self.conflicting_hgnc_ids):
            self.conflicting_hgnc_ids.append( conflicting_hgnc_id )            

        if (not conflicting_hgnc_symbol in self.conflicting_hgnc_symbols):
            self.conflicting_hgnc_symbols.append( conflicting_hgnc_symbol )

        conflicting_db_genes = conflicting_isoform_group.getContainedDatabaseLoci()
        for db, gene_set in conflicting_db_genes.items():
            self.conflicting_db_genes[db].add( gene_set )

        
    def getStrand(self):
        return self.strand
    

    def getHGNCIDAndSymbol(self):
        return (self.hgnc_id, self.hgnc_symbol)


    def getHGNCID(self):
        return self.hgnc_id


    def getHGNCSymbol(self):
        return self.hgnc_symbol


    def getContainedDatabaseLoci(self):
        return self.db_genes


    #def getAltSymbols(self):
    #    return self.gene_symbols


    def hasIsoform(self, isoform):
        return isoform in self.isoforms


    def getIsoforms(self, not_absorbed=False):
        if (not_absorbed):
            return filter(lambda x: not x.wasAbsorbed(), self.isoforms)
        else:
            return self.isoforms


    def getFivePrimeBoundary(self):
        return self.start

    
    def containsIsoformFromSameDBAndLocusAs(self, isoform):
        """Note that UCSD loci are all 'None'"""
        isoform_db, isoform_locus = isoform.getDatabaseAndLocus()
        return (isoform_locus != None and self.db_genes.has_key(isoform_db) and self.db_genes[isoform_db] == isoform_locus)


    def hasDBAndLocusConflict(self, isoform):
        isoform_db, isoform_locus = isoform.getDatabaseAndLocus()
        assert (locus != None)
        return (self.db_genes.has_key(isoform_db) and self.db_genes[isoform_db] != isoform_locus)


    def setGroupNumber(self, group_number):
        self.group_number = group_number


    def setIsoformTSSAndSpliceVariantPolyASiteNumbers(self):
        # Group isoforms by start positions and splice junction usage
        tss_groups = defaultdict(list)
        splice_variant_groups = defaultdict(list)
        polyA_groups = defaultdict(list)
        for isoform in self.getIsoforms(not_absorbed=True): # filter(lambda x: not x.wasAbsorbed(), self.isoforms):
            tss_pos = isoform.getStrandwiseStart()
            tss_groups[tss_pos].append(isoform)

            sj_pos = isoform.getSpliceJunctionPositions()
            sj_pos_tuple = tuple(sorted(sj_pos, reverse=(self.strand == '-')))
            splice_variant_groups[sj_pos_tuple].append(isoform)
            
            polyA_pos = isoform.getStrandwiseStop()
            polyA_groups[polyA_pos].append(isoform)
            
        assert (self.strand != None)

        # TSS numbering increases 5' -> 3' in the -/+ strand context
        tss_positions = tss_groups.keys()
        tss_positions = sorted(tss_positions, reverse=(self.strand == '-'))

        # Set numbering of splice junction patterns
        isoform_splice_variant_number = {}
        for sv_counter, sj_pos_tuple in enumerate(splice_variant_groups.keys()):
            for isoform in splice_variant_groups[sj_pos_tuple]:
                    isoform_splice_variant_number[isoform] = sv_counter + 1
            
        # Set PolyA site numbering, which increases 5' -> 3' in the -/+ strand context
        polyA_positions = polyA_groups.keys()
        polyA_positions = sorted(polyA_positions, reverse=(self.strand == '-'))
        isoform_polyA_variant_number = {}
        for polyA_counter, polyA_pos in enumerate(polyA_positions):
            for isoform in polyA_groups[polyA_pos]:
                isoform_polyA_variant_number[isoform] = polyA_counter + 1
        
        for tss_counter, tss_pos in enumerate(tss_positions):
            for isoform in tss_groups[tss_pos]:
                isoform.setTSSAndSpliceVariantAndPolyASiteNumbers(tss_counter+1, isoform_splice_variant_number[isoform],
                                                                  isoform_polyA_variant_number[isoform])

                
    #def mergeSingleExonIsoforms(self, G, all_isoform_models):
    #    for subG in nx.connected_components(G):
    #        termini = set()
    #        chromosomes = set()
    #        strands = set()
    #        alt_names = []
    #        for isoform in subG:
    #            some_termini = isoform.getExonTerminiPositions()
    #            termini.update(some_termini)
    #            chromosomes.add( isoform.getChromosome() )
    #            strands.add( isoform.getStrand() )
    #            alt_names.append( (isoform, isoform.getSrcDB(), isoform.getSrcDBLocus(), isoform.getTranscriptIDInSrcDB(), isoform.getHGNCID(), isoform.getHGNCSymbol()) )
    #
    #       assert (len(chromosomes) == 1)
    #        chromosome = list(chromosomes)[0]
    #
    #            assert (len(strands) == 1)
    #            strand = list(strands)[0]
    #            
    #            new_start, new_stop = min(termini), max(termini)
    #
    #            new_isoform = RNAIsoform("CGDB", None, None, chromosome, strand)
    #            new_isoform.addExon(chromosome, new_start, new_stop, strand)
    #            new_isoform.setIsoformCompleteness()
    #            for data_tup in alt_names:
    #                new_isoform.addAltNamesData(data_tup)
    #            new_isoform.setLengthAndExonTermini()
    #            new_isoform.setUniqueSignature()
    #
    #            for other_isoform in subG:
    #                new_isoform.updateAbsorptionTracking(other_isoform)  # TODO <-------- Confirm this change 
    #                #new_isoform.recordAbsorbed(other_isoform)
    #            #other_isoform.recordAbsorbedBy(new_isoform)
    #
    #        self.isoforms.add(new_isoform)
    #        all_isoform_models.add(new_isoform)


    #def absorbIsoforms(self, DG):
    #    for subDG in nx.weakly_connected_component_subgraphs(DG):
    #        sorted_subDG = nx.topological_sort(subDG)
    #        root_isoform = sorted_subDG[0]
    #        assert (DG.out_degree(root_isoform) == 0)
    #        for other_isoform in sorted_subDG[1:]:
    #            root_isoform.absorbIsoform(other_isoform)
    

    def performFinalChecks(self, all_isoform_models):
        for isoform in self.isoforms:
            assert (isoform in all_isoform_models)
            isoform.verifyExonTermini()

        # Confirm that there is only one unabsorbed isoform with any signature that is = all exon termini except 3' stop
        termini_signatures = set()
        for isoform in self.getIsoforms(not_absorbed=True):
            termini = isoform.getExonTerminiPositions()
            termini = sorted(termini, reverse=(isoform.getStrand()=='-'))
            sig_tup = tuple(termini[0:-1])

            try:
                assert (sig_tup not in termini_signatures), "IsoformGroup has isoforms with duplicate structures"
            except AssertionError:
                pdb.set_trace()
                
            termini_signatures.add( sig_tup )

        # TODO?: Confirm for every absorbed isoform that its signature (as in 1st TODO) is a subset of some unabsorbed isoform's signature


    def consolidateIsoforms(self, genome_seq): # , all_isoform_models
        # Same value as used in createCGDB.py:standardizePromoterRegion()
        max_allowed_new_orf_len = 24

        isoform_IDs = set(map(lambda x: x.getTranscriptIDInSrcDB(), self.isoforms))

        debug_IDs = set(['CEACAM1.lAug10','ENST00000599389', 'HTR001781.19.1441.8', 'uc010eii', 'uc060zhx'])
        
        # Merge identical isoforms
        G = nx.Graph()
        assert (not any(map(lambda x: x.wasAbsorbed(), self.isoforms))), "ERROR: should not be trying to consolidate already-absorbed isoforms"
        isoforms = list(self.isoforms)
        for i, isoform1 in enumerate(isoforms):
            for isoform2 in isoforms[i+1:]:
                if (set([isoform1.ID,isoform2.ID]) == debug_IDs):
                    pdb.set_trace()
                assert (isoform1 != isoform2)
                if (isoform1.hasSameUniqueSignatureAs(isoform2)):
                    G.add_edge(isoform1, isoform2)

        for identical_isoforms in filter(lambda y: len(y) > 1, nx.clique.find_cliques(G)):
            sorted_identical_isoforms = sorted(identical_isoforms, key=methodcaller("getTranscriptIDInSrcDB")) # For reproducibility
            for other_isoform in sorted_identical_isoforms[1:]:
                sorted_identical_isoforms[0].absorbIsoform(other_isoform, True)


        # Build directed graph of which isoforms can be absorbed by which other isoforms.
        # Graph should be acyclic because identical isoforms already merged.
        unsorted_unabsorbed_isoforms = filter(lambda x: not x.wasAbsorbed(), self.isoforms)
        unabsorbed_isoforms = sorted(unsorted_unabsorbed_isoforms, key=methodcaller("getTranscriptIDInSrcDB")) # For reproducibility
        DG = nx.DiGraph()
        for i, isoform1 in enumerate(unabsorbed_isoforms):
            for isoform2 in self.isoforms:
                if (set([isoform1.ID,isoform2.ID]) <= debug_IDs):
                    print >> sys.stderr, "%s and %s" % (isoform1.ID,isoform2.ID)
                    pdb.set_trace()
                #if (isoform1.getStrandwiseStart() in [141358015,141358016,141358017] and isoform2.getStrandwiseStart() in [141358015,141358016,141358017]):
                #    pdb.set_trace()
                    
                if (isoform1 == isoform2 or G.has_edge(isoform2, isoform1)):
                    continue
                    
                assert (isoform1.getStrand() == isoform2.getStrand())
                isoform1_num_exons = isoform1.getNumberOfExons()
                isoform2_num_exons = isoform2.getNumberOfExons()

                can_merge_1_into_2 = False
                can_merge_2_into_1 = False

                if ((isoform1_num_exons > 1 and isoform2_num_exons == 1) or (isoform1_num_exons == 1 and isoform2_num_exons > 1)):
                    # Absorb if a single-exon isoform is completely contained within exon of the other
                    poten_absorber, poten_absorbee = (isoform1, isoform2) if (isoform1_num_exons > 1) else (isoform2, isoform1)
                    absorbee_start, absorbee_stop = poten_absorbee.getStartStop()
                    
                    for chrom, exon_start, exon_stop in poten_absorber.getExonsAsTuples():
                        if (absorbee_start >= exon_start and absorbee_stop <= exon_stop):
                            if (isoform1_num_exons > 1):
                                can_merge_2_into_1 = True
                            else:
                                can_merge_1_into_2 = True
                            break
                
                elif (isoform1_num_exons > 1 and isoform2_num_exons > 1):
                    isoform1_is_5p_complete = isoform1.isFivePrimeComplete()
                    isoform2_is_5p_complete = isoform2.isFivePrimeComplete()
                    both_5p_complete = isoform1_is_5p_complete == isoform2_is_5p_complete == True

                    isoform1_is_3p_complete = isoform1.isThreePrimeComplete()
                    isoform2_is_3p_complete = isoform2.isThreePrimeComplete()
                    both_3p_complete = isoform1_is_3p_complete == isoform2_is_3p_complete == True
                
                    isoform1_5p_start = isoform1.getStrandwiseStart()
                    isoform2_5p_start = isoform2.getStrandwiseStart()
                    
                    isoform1_3p_stop = isoform1.getStrandwiseStop()
                    isoform2_3p_stop = isoform2.getStrandwiseStop()

                    same_5p_pos = isoform1_5p_start == isoform2_5p_start
                    same_3p_pos = isoform1_3p_stop == isoform2_3p_stop

                    isoform1_SJs = isoform1.getIntronTuples()
                    isoform2_SJs = isoform2.getIntronTuples()

                    isoform1_exons = isoform1.getExonTuples()
                    isoform2_exons = isoform2.getExonTuples()
                    
                    #if (len(set([isoform1_5p_start,isoform2_5p_start]) & set([141367265,141367267]))>0):
                    #    pdb.set_trace()

                    # The following logic allows mRNAs with incomplete ends to be absorbed.
                    if (isoform1_SJs == isoform2_SJs):
                        can_merge_1_into_2_SJ = True
                        can_merge_2_into_1_SJ = True

                        if (same_5p_pos):
                            can_merge_1_into_2_5p = True
                            can_merge_2_into_1_5p = True
                        else:
                            if (isoform2.firstExonOverlaps(isoform1_5p_start)):
                                can_merge_1_into_2_5p = isoform1_is_5p_complete!=True or not isoform1.fivePrimeExtensionCreatesORF(isoform2_5p_start, genome_seq, max_allowed_new_orf_len)
                                can_merge_2_into_1_5p = False #isoform2_is_5p_complete!=False and not isoform2.fivePrimeContractionRemovesORF(isoform1_5p_start, genome_seq)
                            else:
                                can_merge_2_into_1_5p = isoform2_is_5p_complete!=True or not isoform2.fivePrimeExtensionCreatesORF(isoform1_5p_start, genome_seq, max_allowed_new_orf_len)
                                can_merge_1_into_2_5p = False #isoform1_is_5p_complete!=False and not isoform1.fivePrimeContractionRemovesORF(isoform2_5p_start, genome_seq)
                                
                        # This logic allows an isoform with a longer 3' UTR to be absorbed by an isoform 
                        # with a shorter 3' UTR if  the shorter has it's polyA site within the longer's 3' UTR.
                        # Make certain that merging isn't bidirectional. If they are at this point,
                        # make the isoform with the *shorter* 3' UTR be the absorbing isoform.
                        can_merge_1_into_2_3p = same_3p_pos or isoform2.lastExonOverlaps(isoform1_3p_stop)
                        can_merge_2_into_1_3p = same_3p_pos or isoform1.lastExonOverlaps(isoform2_3p_stop)

                        #if (can_merge_1_into_2_5p and can_merge_2_into_1_5p and can_merge_1_into_2_3p and can_merge_2_into_1_3p):
                        #    if (same_3p_pos):
                        #        # Use the isoform with the shorter first exon as the absorber
                        #        (can_merge_2_into_1_3p, can_merge_1_into_2_3p) = (False, True) if (isoform1.getFirstExon().getLength() > isoform2.getFirstExon().getLength()) else (True, False)
                        #    else:
                        #        # Use the isoform with the shorter last exon as the absorber
                        #        (can_merge_2_into_1_3p, can_merge_1_into_2_3p) = (False, True) if (isoform1.getLastExon().getLength() > isoform2.getLastExon().getLength()) else (True, False)

                    elif (isoform1_SJs.issubset(isoform2_SJs)):
                        # Overlap if (start1 <= stop2 and start2 <= stop1)
                        # True if no intron of isoform2 overlaps an exon of isoform 1
                        can_merge_1_into_2_SJ = not( any( map(lambda x: x[0][0]<=x[1][1] and x[1][0]<=x[0][1], list(product(isoform2_SJs,isoform1_exons))) ) )

                        if (can_merge_1_into_2_SJ and isoform2.doesOverlapPosition(isoform1_5p_start) and isoform2.doesOverlapPosition(isoform1_3p_stop)):
                            # If isoform1 is completely contained within isoform2, then absorb isoform1 -> isoform2
                            can_merge_1_into_2_5p = True
                            can_merge_1_into_2_3p = True
                        elif (can_merge_1_into_2_SJ and not isoform2.firstExonOverlaps(isoform1_5p_start) and not isoform2.lastExonOverlaps(isoform1_3p_stop) and
                            isoform1_is_5p_complete!=True and isoform1_is_3p_complete!=True):
                            can_merge_1_into_2_5p = True
                            can_merge_1_into_2_3p = True
                        else:
                            if (same_5p_pos or (isoform2.doesOverlapPosition(isoform1_5p_start) and isoform1_is_5p_complete==False)):
                                can_merge_1_into_2_5p = True
                            elif (isoform2.firstExonOverlaps(isoform1_5p_start)):
                                can_merge_1_into_2_5p = not isoform1.fivePrimeExtensionCreatesORF(isoform2_5p_start, genome_seq, max_allowed_new_orf_len)
                            #elif (isoform1.firstExonOverlaps(isoform2_5p_start)):
                            #    can_merge_1_into_2_5p = not isoform1.fivePrimeContractionRemovesORF(isoform2_5p_start, genome_seq)
                            else:
                                can_merge_1_into_2_5p = False

                            # This logic allows an isoform with a longer 3' UTR to be absorbed by an isoform with a shorter 3' UTR if the shorter has it's
                            # polyA site within the longer's 3' UTR and longer is 3' complete (ie doesn't extend to unknown exons)
                            can_merge_1_into_2_3p = (isoform2.doesOverlapPosition(isoform1_3p_stop) and isoform1_is_3p_complete!=True) #or \
                                                    #((isoform1.lastExonOverlaps(isoform2_3p_stop) or isoform2.lastExonOverlaps(isoform1_3p_stop)) and both_3p_complete)

                        can_merge_2_into_1_SJ = False
                        can_merge_2_into_1_5p = False # Doesn't matter, since can_merge_2_into_1_SJ is False
                        can_merge_2_into_1_3p = False
                        
                    elif (isoform2_SJs.issubset(isoform1_SJs)):
                        # True if no intron of isoform1 overlaps an exon of isoform 2
                        can_merge_2_into_1_SJ = not( any( map(lambda x: x[0][0]<=x[1][1] and x[1][0]<=x[0][1], list(product(isoform1_SJs,isoform2_exons))) ) )

                        if (can_merge_2_into_1_SJ and isoform1.doesOverlapPosition(isoform2_5p_start) and isoform1.doesOverlapPosition(isoform2_3p_stop)):
                            # If isoform2 is completely contained within isoform1, then absorb isoform2 -> isoform1
                            can_merge_2_into_1_5p = True
                            can_merge_2_into_1_3p = True
                        elif (can_merge_2_into_1_SJ and not isoform1.firstExonOverlaps(isoform2_5p_start) and not isoform1.lastExonOverlaps(isoform2_3p_stop) and
                            isoform2_is_5p_complete!=True and isoform2_is_3p_complete!=True):
                            can_merge_2_into_1_5p = True
                            can_merge_2_into_1_3p = True
                        else:
                            if (same_5p_pos or (isoform1.doesOverlapPosition(isoform2_5p_start) and isoform2_is_5p_complete==False)):
                                can_merge_2_into_1_5p = True
                            elif (isoform1.firstExonOverlaps(isoform2_5p_start)):
                                can_merge_2_into_1_5p = not isoform2.fivePrimeExtensionCreatesORF(isoform1_5p_start, genome_seq, max_allowed_new_orf_len)
                            #elif (isoform2.firstExonOverlaps(isoform1_5p_start)):
                            #    can_merge_2_into_1_5p = not isoform2.fivePrimeContractionRemovesORF(isoform1_5p_start, genome_seq)
                            else:
                                can_merge_2_into_1_5p = False
                            
                            can_merge_2_into_1_3p = (isoform1.doesOverlapPosition(isoform2_3p_stop) and isoform2_is_3p_complete!=True) #or \
                                                    #((isoform1.lastExonOverlaps(isoform2_3p_stop) or isoform2.lastExonOverlaps(isoform1_3p_stop)) and both_3p_complete)

                        can_merge_1_into_2_SJ = False
                        can_merge_1_into_2_5p = False
                        can_merge_1_into_2_3p = False

                    else:
                        can_merge_1_into_2_5p = False
                        can_merge_2_into_1_5p = False
                        can_merge_1_into_2_SJ = False
                        can_merge_2_into_1_SJ = False
                        can_merge_1_into_2_3p = False
                        can_merge_2_into_1_3p = False

                    can_merge_1_into_2 = can_merge_1_into_2_5p and can_merge_1_into_2_3p and can_merge_1_into_2_SJ
                    can_merge_2_into_1 = can_merge_2_into_1_5p and can_merge_2_into_1_3p and can_merge_2_into_1_SJ

                elif (isoform1_num_exons == 1 and isoform2_num_exons == 1):
                    start1, stop1 = isoform1.getStartStop()
                    start2, stop2 = isoform2.getStartStop()

                    can_merge_1_into_2 = start1>=start2 and stop1<=stop2
                    can_merge_2_into_1 = start2>=start1 and stop2<=stop1


                assert (not(can_merge_1_into_2 and can_merge_2_into_1)), "ERROR: should not be able to merge both isoforms into each other at this point"
                    
                if (can_merge_1_into_2 or can_merge_2_into_1):
                    absorber, absorbee = (isoform2, isoform1) if (can_merge_1_into_2) else (isoform1, isoform2)
                    DG.add_edge(absorbee, absorber) # Edge (isoformA, isoformB) means absorb isoformB into isoformB. Or, isoformA -> isoformB

        # Absorb into each node with out_degree==0 all isoforms that are reachable from it
        assert (nx.is_directed_acyclic_graph(DG)), "ERROR: absorption graph is not acyclic"
        for absorber_isoform in filter(lambda x: DG.out_degree(x) == 0, DG.nodes()):
            for absorbee_isoform in nx.ancestors(DG, absorber_isoform):
                assert (absorber_isoform != absorbee_isoform)
                if (not absorber_isoform.hasAlreadyAbsorbed(absorbee_isoform)):
                    absorber_isoform.absorbIsoform(absorbee_isoform, False)


        # Merge overlapping single-exon isoforms
        #unabsorbed_isoforms = filter(lambda x: not x.wasAbsorbed(), self.isoforms)
        #G_single_exon_merge = nx.Graph() # Merge the isoforms in each connected component
        #for i, isoform1 in enumerate(unabsorbed_isoforms):
        #    for isoform2 in unabsorbed_isoforms[i+1:]:
        #        if (isoform1.getNumberOfExons()==1 and isoform2.getNumberOfExons()==1):
        #            start1, stop1 = isoform1.getStartStop()
        #            start2, stop2 = isoform2.getStartStop()
        #            if (start1 <= stop2 and start2 <= stop1): # Do overlap
        #                G_single_exon_merge.add_edge(isoform1, isoform2)
        #self.mergeSingleExonIsoforms(G_single_exon_merge, all_isoform_models)

        
    def setLocusSymbol(self):
        self.locus_symbol = "None"
        if (self.hgnc_symbol != None):
            self.locus_symbol = self.hgnc_symbol
        else:
            for db in ["RefSeq", "GENCODE", "HInv", "AceView"]:
                locus_symbols_from_isoforms = set()
                for isoform in self.isoforms:
                    isoform_locus_symbol = isoform.getLocusSymbolFromDatabase(db)
                    if (isoform_locus_symbol != None):
                        locus_symbols_from_isoforms.add( isoform_locus_symbol )
                if (len(locus_symbols_from_isoforms) == 1):
                    self.locus_symbol = list(locus_symbols_from_isoforms)[0]
                    break
                elif (len(locus_symbols_from_isoforms) > 1):
                    self.locus_symbol = ",".join(list(locus_symbols_from_isoforms))
                    print >> sys.stderr, "WARNING: in naming isoform group, multiple %s locus symbols from isoforms: %s." % \
                        (db, " ".join(locus_symbols_from_isoforms))

                    
    def writeAsBED(self, op, write_absorbed):
        color_default = "0,0,0"
        color_for_absorbed = "255,0,0"

        for isoform in self.isoforms:
            if (not isoform.wasAbsorbed()):
                bed_fields = isoform.getAsBED12Fields("CGDB", None, color_default)
                op.write("%s\n" % "\t".join(bed_fields))
            elif (write_absorbed):
                bed_fields = isoform.getAsBED12Fields("ID", None, color_for_absorbed)
                op.write("%s\n" % "\t".join(bed_fields))


    def writeAsGTF(self, op, write_absorbed):
        score_default = 1000
        score_for_absorbed = 500

        for isoform in self.isoforms:
            if (not isoform.wasAbsorbed()):
                gtf_fields = isoform.getAsGTF("CGDB", score_default)
                op.write("%s\n" % "\n".join(gtf_fields))
            elif (write_absorbed):
                gtf_fields = isoform.getAsGTF("ID", score_for_absorbed)
                op.write("%s\n" % "\n".join(gtf_fields))


    def writeAsGFF3(self, op, write_absorbed):
        ig_gff3_ID = "ig%d" % self.group_number

        op.write("%s\tCGDB\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=CGDB%d;gene=%s\n" % \
                 (self.chromosome, self.start, self.stop, self.strand, ig_gff3_ID, self.group_number, self.locus_symbol))

        all_mRNA_tuples = []
        exon_parents = defaultdict(list)
        for isoform in self.getIsoforms(not_absorbed=not write_absorbed):
            exon_parent = isoform.getCGDBName()
            mRNA_tuple, exon_tuples = isoform.getAsGFF3(ig_gff3_ID)
            all_mRNA_tuples.append(mRNA_tuple)
            for exon_tuple in exon_tuples:
                exon_parents[exon_tuple].append(exon_parent)
            
        all_mRNA_tuples = sorted(all_mRNA_tuples, key=itemgetter(1), reverse=(self.strand == '-'))
        for mRNA_tuple in all_mRNA_tuples:
            op.write("%s\tCGDB\tmRNA\t%d\t%d\t.\t%s\t.\t%s\n" % mRNA_tuple)
            
        all_exon_tuples = sorted(exon_parents.keys(), key=itemgetter(1), reverse=(self.strand == '-'))
        for exon_tuple in all_exon_tuples:
            parents = ",".join(exon_parents[exon_tuple])
            new_exon_tuple = tuple(list(exon_tuple) + ["Parent=%s;" % parents])
            op.write("%s\tCGDB\texon\t%d\t%d\t.\t%s\t.\t%s\n" % new_exon_tuple)

        op.write("###\n")
