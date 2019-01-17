from __future__ import print_function

from collections import defaultdict
from itertools import combinations
import numpy as np
from operator import attrgetter, methodcaller

from Bio.Seq import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
my_codon_table = CodonTable.unambiguous_dna_by_id[1]

# Scientific consensus is that the full set of alternative start codons in vivo is ATG, CTG, TTG, GTG, ACG, ATT, ATC, ATA, AGG, AAG
#my_codon_table.start_codons.extend( ["ACG", "GTG"] ) <- UNCOMMENT WHEN UPDATED

try:
    DNA_complement_table = str.maketrans("ACGTNacgtn","TGCANtgcan")
except AttributeError:  # TODO: Remove when all upgraded to python3
    import string
    DNA_complement_table = string.maketrans("ACGTNacgtn","TGCANtgcan")

import sys
from pyfaidx import Fasta
import pdb


class Exon(object):
    def __init__(self, parent_isoform_instance, chromosome, start, stop, strand):
        self.parent_isoform = parent_isoform_instance
        self.chromosome = chromosome
        self.strand = strand

        # These two are relative to the '+' strand, irrespective of the strand on which the exon actually lies
        self.start = start
        self.stop = stop

        self.orig_start = start
        self.orig_stop = stop

        # These are only set for 1st exons. They indicate the most downstream position from the
        # strandwise exon start for where there is evidence of the start of transcription
        self.original_five_prime_position = None
        self.most_downstream_robust_start = None
        self.most_downstream_permissive_start = None
        
        self.is_five_prime_complete = None
        self.is_three_prime_complete = None

        
    def switchStrand(self):
        self.is_five_prime_complete = None
        self.is_three_prime_complete = None
        self.strand = '+' if (self.strand == '-') else '-'


    def getLengthAdjustmentAmounts(self):
        if (self.strand == '+'):
            length_adjustment_5p = self.orig_start - self.start
            length_adjustment_3p = self.stop - self.orig_stop
        else:
            length_adjustment_5p = self.stop - self.orig_stop
            length_adjustment_3p = self.orig_start - self.start

        return (length_adjustment_5p, length_adjustment_3p)


    def getSrcDB(self):
        return self.parent_isoform.getSrcDB()


    def getStrand(self):
        return self.strand


    def getStart(self):
        return self.start


    def getStop(self):
        return self.start


    def getStartStop(self):
        return (self.start, self.stop)
    

    def getLength(self):
        return self.stop - self.start + 1


    def getStrandwiseStart(self):
        return self.start if (self.strand == '+') else self.stop


    def getStrandwiseStop(self):
        return self.stop if (self.strand == '+') else self.start


    def getOriginalStrandwiseStart(self):
        assert (self.original_five_prime_position != None) # This is only set for 1st exons
        return self.original_five_prime_position
        

    def setAsFirstExon(self):
        if (self.strand == '-'):
            self.original_five_prime_position = self.stop
            self.most_downstream_robust_start = self.stop
            self.most_downstream_permissive_start = self.stop
        else:
            self.original_five_prime_position = self.start
            self.most_downstream_robust_start = self.start
            self.most_downstream_permissive_start = self.start


    def unsetAsFirstExon(self):
        if (self.strand == '-'):
            self.original_five_prime_position = None
            self.most_downstream_robust_start = None
            self.most_downstream_permissive_start = None
        else:
            self.original_five_prime_position = None
            self.most_downstream_robust_start = None
            self.most_downstream_permissive_start = None


    def isFirstExon(self):
        return (self.original_five_prime_position != None)
    

    def setStrandwiseStart(self, pos):
        if (self.strand == '-'):
            self.stop = pos
        else:
            self.start = pos
            

    def setStrandwiseStop(self, pos):
        if (self.strand == '+'):
            self.stop = pos
        else:
            self.start = pos


    def resetStrandwiseStart(self, pos):
        try:
            if (self.strand == '-'):
                assert (pos > self.start or self.start==self.stop), "Strandwise start for exon being reset outside (downstream) of exon"
                self.setIfMostDownstreamVariableStart(min(pos,self.stop), "robust")
                self.stop = pos
            else:
                assert (pos < self.stop or self.start==self.stop), "Strandwise start for exon being reset outside (downstream) of exon"
                self.setIfMostDownstreamVariableStart(max(pos,self.start), "robust")
                self.start = pos
        except AssertionError as ae:
            pdb.set_trace()
            print(ae, file=sys.stderr)


    def resetStrandwiseStop(self, pos):
        if (self.strand == '+'):
            assert (pos > self.start), "Strandwise stop for exon being reset outside (upstream) of exon"
            self.stop = pos
        else:
            assert (pos < self.stop), "Strandwise stop for exon being reset outside (upstream) of exon"
            self.start = pos


    def setIfMostDownstreamVariableStart(self, pos, label):
        assert (self.original_five_prime_position != None), "Setting variable start on exon that is not a first exon"
        assert (pos >= self.start and pos <= self.stop), "New variable start position on first exon is not within exon"
        
        if (label == "robust"):
            if ((self.strand == '-' and pos < self.most_downstream_robust_start) or
                (self.strand == '+' and pos > self.most_downstream_robust_start)):
                self.most_downstream_robust_start = pos
        elif (label == "permissive"):
            if ((self.strand == '-' and pos < self.most_downstream_permissive_start) or
                (self.strand == '+' and pos > self.most_downstream_permissive_start)):
                self.most_downstream_permissive_start = pos
        else:
            print("ERROR: unrecognized variable start label -> %s. Exiting." % label, file=sys.stderr)
            sys.exit(1)
            
    def getStop(self):
        return self.stop
    
    def getStartStop(self):
        return (self.start, self.stop)
    
    def getFivePrimePosition(self):
        return self.start if (self.strand == '+') else self.stop

    def getThreePrimePosition(self):
        return self.stop if (self.strand == '+') else self.start

    def getStrandwiseSequence(self, genome_seq, as_string=False):
        """Returns the nucleotide sequence (as a Seq object) of the exon as it would be read during transcription.
        That is, for '-' strand exon the nucleotide sequence returned is the reverse complement
        of the genome reference sequence."""
        #region_spec = "%s:%d-%d" % (self.chromosome, self.start, self.stop)
        nuc_seq_fasta = genome_seq[self.chromosome][self.start-1:self.stop] #pysam.faidx(genome_seq, region_spec)
        nuc_seq = str(nuc_seq_fasta) # ''.join(map(lambda x: x.strip(), nuc_seq_fasta[1:]))
        #nuc_seq = nuc_seq.upper()
        #s = Seq(nuc_seq, IUPAC.unambiguous_dna)
        if (self.strand == '-'):
            nuc_seq = nuc_seq[::-1].translate(DNA_complement_table)
            
        #if (as_string):
        #    s = str(s)
    
        return nuc_seq


    def setFivePrimeIncomplete(self):
        self.is_five_prime_complete = False


    def setThreePrimeIncomplete(self):
        self.is_three_prime_complete = False


    def setFivePrimeComplete(self):
        self.is_five_prime_complete = True


    def setThreePrimeComplete(self):
        self.is_three_prime_complete = True


    def isFivePrimeComplete(self):
        return self.is_five_prime_complete


    def isThreePrimeComplete(self):
        return self.is_three_prime_complete


    def isDownstreamOf(self, other_exon):
        assert self.strand == other_exon.getStrand()
        return ((self.strand == '+' and self.start > other_exon.getStop()) or
                (self.strand == '-' and self.stop < other_exon.getStart()))


    def doesOverlapExon(self, other_exon):
        assert self.strand == other_exon.getStrand()
        return not(other_exon.getStop() < self.start or other_exon.getStart() > self.stop)


    def doesOverlapPosition(self, pos):
        return not(pos < self.start or pos > self.stop)

        
    def fivePrimeContractionRemovesORF(self, downstream_start, genome_seq):
        assert (self.doesOverlapPosition(downstream_start)), "Trying to contract exon to a 5' position not with the exon"
        if (self.strand == '+'):
            remove_nuc_seq = genome_seq[self.chromosome][self.start:downstream_start+1].seq.upper()
        else:
            remove_nuc_seq = genome_seq[self.chromosome][downstream_start:self.stop].seq.upper()
            remove_nuc_seq = str(remove_nuc_seq)
            try:
                remove_nuc_seq = remove_nuc_seq[::-1].translate(DNA_complement_table)
            except TypeError as te:
                pdb.set_trace()
                
        return any(map(lambda x: x in remove_nuc_seq, my_codon_table.start_codons))
        

    def fivePrimeExtensionCreatesORF(self, upstream_start, genome_seq, max_allowed_new_orf_len):
        if (self.strand == '+'):
            exten_nuc_seq = genome_seq[self.chromosome][upstream_start-1: self.start-1].seq.upper()
        else:
            exten_nuc_seq = genome_seq[self.chromosome][self.stop:upstream_start].seq.upper()

        exten_nuc_seq = Seq(exten_nuc_seq, IUPAC.unambiguous_dna)
        if (self.strand == '-'):
            exten_nuc_seq = exten_nuc_seq.reverse_complement()

        isoform_nuc_seq = None
        new_orf_is_created = False
        
        # Determine whether any potentially significant ORFs are initiated in the extension region
        for frame in [0,1,2]:
            complete_codons = exten_nuc_seq[frame:len(exten_nuc_seq)-(len(exten_nuc_seq[frame:])%3)]
            trailing_partial_codon = exten_nuc_seq[len(exten_nuc_seq)-(len(exten_nuc_seq[frame:])%3):]
            
            try:
                translations = complete_codons.translate(my_codon_table)
            except TranslationError as te:
                if ('N' in complete_codons or 'n' in complete_codons):
                    new_orf_is_created = True
                    break
                else:
                    import pdb
                    pdb.set_trace()
                
            orfs = translations.split('*')
            if (complete_codons != '' and len(orfs) > 0):
                if ('M' in orfs[-1]):
                    exten_orf_len = len(orfs[-1]) - orfs[-1].find('M')
                    if (exten_orf_len >= max_allowed_new_orf_len):
                        new_orf_is_created = True
                    elif (translations[-1] != '*'):
                        # Determine if a long ORF that extends into original mRNA is created
                        if (isoform_nuc_seq == None):
                            isoform_nuc_seq = self.parent_isoform.getStrandwiseSequence(3*max_allowed_new_orf_len, genome_seq)

                        composite_nuc_seq = trailing_partial_codon  + isoform_nuc_seq
                        composite_orf = composite_nuc_seq[0:len(composite_nuc_seq)-(len(composite_nuc_seq)%3)]
                        composite_orf_translation = composite_orf.translate(my_codon_table)
                        if (exten_orf_len + len(composite_orf_translation.split('*')[0]) >= max_allowed_new_orf_len):
                            new_orf_is_created = True
                            break
                        
                if (len(orfs) > 1):
                    for orf in filter(lambda x: 'M' in x, orfs[0:-1]):
                        orf_len = len(orf) - orf.find('M')
                        if (orf_len >= max_allowed_new_orf_len):
                            new_orf_is_created = True
                            break

            if (new_orf_is_created):
                break

        return new_orf_is_created


    def getAsGFF3(self):
        return (self.chromosome, self.start, self.stop, self.strand)


class RNAIsoform(object):
    def __init__(self, src_db, locus_ID, transcript_ID, chromosome, strand):
        self.tss_number = None
        self.splice_variant_number = None
        self.polyA_variant_number = None
        self.CGDB_name = None
        
        self.src_db = src_db
        self.chromosome = chromosome
        self.src_db_locus = locus_ID
        self.strand = strand
        self.ID = transcript_ID
        self.length = None
        self.nuc_seq = None

        self.is_five_prime_complete = None
        self.is_three_prime_complete = None
        
        self.hgnc_id = None
        self.hgnc_symbol = None

        self.alt_names = []
        self.polyA_sites = set() # These should all be internal to the last exon
        
        # Format: self.locus_names[database] -> name and synonyms obtained from database, in preference order
        self.locus_names = {}

        # To record with isoform(s) absorbed this isoform, and which isoform(s) this isoform absorbed.
        self.absorbed_by = set()
        self.absorbed = set()
        
        # Exons are in the order first-to-last if isoform is on '+' strand,
        # and in the order last-to-first if on the '-' strand
        self.ordered_exons = []
        self.exon_termini_positions = None
        self.unique_signature = None

        self.coords_mRNA_to_genome = None
        self.coords_genome_to_mRNA = None
        self.genomic_positions = None

        #self.is_known_fusion = False


    def __lt__(self, other):
        return self.getCGDBName() < other.getCGDBName()


#    def setAsFusionTranscript(self):
#        self.is_known_fusion = True


#    def isKnownFusionTranscript(self):
#        return self.is_known_fusion


    def fivePrimeContractionRemovesORF(self, downstream_start, genome_seq):
        first_exon = self.getFirstExon()
        return first_exon.fivePrimeContractionRemovesORF(downstream_start, genome_seq)


    def fivePrimeExtensionCreatesORF(self, upstream_start, genome_seq, max_allowed_new_orf_len):
        first_exon = self.getFirstExon()
        return first_exon.fivePrimeExtensionCreatesORF(upstream_start, genome_seq, max_allowed_new_orf_len)


    def switchStrand(self):
        self.strand = '+' if (self.strand == '-') else '-'
        self.is_five_prime_complete, self.is_three_prime_complete = (self.is_three_prime_complete, self.is_five_prime_complete)

        first_exon = self.getFirstExon()
        first_exon.unsetAsFirstExon()

        for exon in self.ordered_exons:
            exon.switchStrand()

        self.polyA_sites.clear()
        self.setIsoformCompleteness()
        
        assert (len(self.absorbed_by)==0)
        assert (len(self.absorbed)==0)
        assert (self.coords_mRNA_to_genome == None)
        assert (self.coords_genome_to_mRNA == None)
        assert (self.genomic_positions == None)
        assert (self.nuc_seq == None)
        
        
    def verifyExonTermini(self):
        first_exon_index = 0 if (self.strand == '+') else -1
        last_exon_index = -1 if (self.strand == '+') else 0

        # Verify Completeness consistency between isoform and exons
        assert (self.is_five_prime_complete == self.ordered_exons[first_exon_index].isFivePrimeComplete())
        assert (self.is_three_prime_complete == self.ordered_exons[last_exon_index].isThreePrimeComplete())

        # Verify all internal exons are 5' and 3' complete
        if (len(self.ordered_exons) > 2):
            for exon in self.ordered_exons[1:-1]:
                try:
                    assert (exon.isFivePrimeComplete() == exon.isThreePrimeComplete() != False)
                except AssertionError:
                    import pdb
                    pdb.set_trace()
                    
        if (len(self.ordered_exons) > 1):
            # Verify first exon is 3' complete and last exon is 5' complete
            assert (self.ordered_exons[first_exon_index].isThreePrimeComplete() != False)
            assert (self.ordered_exons[last_exon_index].isFivePrimeComplete() != False)
            
        # Verify that 3' position of last exon is not a polyA site if exon is 3' incomplete
        if (self.is_three_prime_complete == False):
            try:
                assert (self.getStrandwiseStop() not in self.polyA_sites), "Incomplete 3' end used as polyA site"
            except AssertionError:
                pdb.set_trace()
                
        # Verify that the termini of the exons are what is actually recorded as exon termini for the isoform
        exon_termini = set()
        has_1bp_exon = False
        for e in self.ordered_exons:
            exon_start, exon_stop = e.getStartStop()
            if (exon_start == exon_stop):
                has_1bp_exon = True  # See, for instance, isoform ENST00000619410.4 of SCN7A
            exon_termini.add(exon_start)
            exon_termini.add(exon_stop)

        try:
            assert (exon_termini == self.exon_termini_positions)
            assert ((len(self.exon_termini_positions) >= 2 and len(self.exon_termini_positions) % 2 == 0) or has_1bp_exon), \
              "ERROR: In verifyExonTermini() - malformed isoform, problem with number of exon termini."
        except AssertionError as ae:
            pdb.set_trace()


    def absorbIsoform(self, other_isoform, is_identical):
        '''Does not delete other_isoform.'''

        if (is_identical):
            assert (self.exon_termini_positions == other_isoform.exon_termini_positions)

        # If one isoform is annotated incomplete and the other is not, then defer to the "complete" isoform. Exception is when the other isoform
        # derives from the AceView database, as its downloadable files don't contain the completeness information and so could inaccurately indicate
        # completeness.
        #if ((not self.isFivePrimeComplete() and other_isoform.isFivePrimeComplete() and other_isoform.getSrcDB() != "AceView") and
        #    (self.getStrandwiseStart() == other_isoform.getStrandwiseStart())):
        #    self.setFivePrimeComplete()
        #    first_exon_index = 0 if (self.strand == '+') else -1
        #    self.ordered_exons[first_exon_index].setFivePrimeComplete()

        #if ((not self.isThreePrimeComplete() and other_isoform.isThreePrimeComplete() and other_isoform.getSrcDB() != "AceView") and
        #    (self.getStrandwiseStop() == other_isoform.getStrandwiseStop())):
        #    self.setThreePrimeComplete()
        #    last_exon_index = -1 if (self.strand == '+') else 0
        #    self.ordered_exons[last_exon_index].setThreePrimeComplete()
        #    stop_pos = self.getStrandwiseStop()
        #    self.addPolyASite(stop_pos)
        #    #other_stop_pos = other_isoform.getStrandwiseStop() 
        #    #other_isoform.discardPolyASite(other_stop_pos)
            
        assert (self.strand == other_isoform.getStrand())
        
        # Merge HGNC ID & Symbol
        if (self.hgnc_id == None and other_isoform.hgnc_id != None):
            self.hgnc_id = other_isoform.hgnc_id
            self.hgnc_symbol = other_isoform.hgnc_symbol
        elif (self.hgnc_id != other_isoform.hgnc_id and other_isoform.hgnc_id != None):
            self.alt_names.append( (other_isoform.src_db, other_isoform.src_db_locus, other_isoform.ID, other_isoform.hgnc_id, other_isoform.hgnc_symbol) )
                
        for db, names in other_isoform.locus_names.items():
            if (db not in self.locus_names):
                self.locus_names[db] = names
            else:
                for name in names:
                    if (not name in self.locus_names[db]):
                        self.locus_names[db].append(name)

        # Add alt polyA site if necessary
        if (not is_identical):
            this_last_exon = self.getLastExon()
            other_last_exon = other_isoform.getLastExon()
            if (this_last_exon.getStrandwiseStart() == other_last_exon.getStrandwiseStart()): # If same last exon...
                this_stop = this_last_exon.getStrandwiseStop()
                other_stop = other_last_exon.getStrandwiseStop()
                if ((self.strand == '+' and this_stop > other_stop) or (self.strand == '-' and this_stop < other_stop)):
                    if (other_isoform.isThreePrimeComplete()):
                        self.addPolyASite(other_stop)  # The polyA site of the other isoforms becomes an alt internal (or identical) polyA site
                elif (this_stop != other_stop):
                    # The absorbee has the longer 3' UTR.
                    self.addPolyASite(other_stop)
                    
        other_isoform.recordAbsorbedBy(self)
        self.recordAbsorbed(other_isoform)


    def isInsideRange(self, chromosome, left_position, right_position):
        return (self.chromosome==chromosome) and (min(self.exon_termini_positions) >= left_position) and (min(self.exon_termini_positions) <= right_position)


    def firstExonOverlaps(self, pos):
        e = self.getFirstExon()
        return e.doesOverlapPosition(pos)


    def lastExonOverlaps(self, pos):
        e = self.getLastExon()
        return e.doesOverlapPosition(pos)

        
    def updateAbsorptionTracking(self, other_isoform):
        # Existing 'absorbed' and 'absorbed_by' information in other_isoform needs to be transferred to this isoform
        other_absorbed = other_isoform.getAbsorbed()
        self.absorbed.update(other_absorbed)
        other_isoform.resetAbsorbed()

        #other_absorbed_by = other_isoform.getAbsorbedBy()
        #self.absorbed_by.update(other_absorbed_by)
        #other_isoform.resetAbsorbedBy()

        #for isoform in other_absorbed_by:
        #    if (isoform != self):
        #        isoform.removeFromAbsorbed(other_isoform)
        #        other_isoform.removeFromAbsorbedBy(isoform)
            
        other_isoform.recordAbsorbedBy(self)
        self.recordAbsorbed(other_isoform)


    def recordAbsorbedBy(self, absorber_isoform):
        assert (not absorber_isoform in self.absorbed_by), "ERROR: this isoform should only be absorbed by a particular isoform once"
        self.absorbed_by.add(absorber_isoform)


    def recordAbsorbed(self, absorbee_isoform):
        assert (not absorbee_isoform in self.absorbed), "ERROR: this isoform should only absorb a particular isoform once"
        self.absorbed.add(absorbee_isoform)


    def hasAlreadyAbsorbed(self, an_isoform):
        return an_isoform in self.absorbed

        
    def removeFromAbsorbed(self, an_isoform):
        self.absorbed.remove(an_isoform)


    def removeFromAbsorbedBy(self, an_isoform):
        self.absorbed_by.remove(an_isoform)


    def resetAbsorbed(self):
        self.absorbed.clear()


    def resetAbsorbedBy(self):
        self.absorbed_by.clear()


    def wasAbsorbed(self):
        return (len(self.absorbed_by) > 0)

        
    def getAbsorbed(self):
        return self.absorbed


    def getAbsorbedBy(self):
        return self.absorbed_by
    

    def addAltNamesData(self, data_tup):
        '''Tuple is of the form (src DB, src DB locus, transcript ID, HGNC ID, HGNC Symbol)'''
        self.alt_names.append(data_tup)

        
    def setIsoformCompleteness(self):
        """For any needed actions once the isoform has been completely read/instantiated."""
        first_exon_index = 0 if (self.strand == '+') else -1
        last_exon_index = -1 if (self.strand == '+') else 0

        # Very short first exons in multi-exon isoforms are most likely 5' incomplete, unless annotated otherwise
        if (len(self.ordered_exons)>1 and self.ordered_exons[first_exon_index].getLength() <= 10 and self.is_five_prime_complete==None):
            self.is_five_prime_complete = False

        if (self.is_five_prime_complete == True):
            self.ordered_exons[first_exon_index].setFivePrimeComplete()
        elif (self.is_five_prime_complete == False):
            self.ordered_exons[first_exon_index].setFivePrimeIncomplete()

        if (self.is_three_prime_complete == True):
            self.ordered_exons[last_exon_index].setThreePrimeComplete()
            polyA_site = self.ordered_exons[last_exon_index].getStrandwiseStop()
            self.addPolyASite( polyA_site )
        elif (self.is_three_prime_complete == False):
            self.ordered_exons[last_exon_index].setThreePrimeIncomplete()

        first_exon = self.getFirstExon()
        first_exon.setAsFirstExon()

        # Check that exons are all separate
        exons_dont_overlap = not any(map(lambda exons: exons[0].doesOverlapExon(exons[1]), combinations(self.ordered_exons,2)))

        return exons_dont_overlap

    
    def setUniqueSignature(self):
        assert (len(self.exon_termini_positions) > 0)
        self.unique_signature = ",".join(([self.chromosome, self.strand] + map(str, sorted(self.exon_termini_positions))))


    def doesOverlapPosition(self, pos):
        return any(map(lambda x: x.doesOverlapPosition(pos), self.ordered_exons))

        
    def getUniqueSignature(self):
        return self.unique_signature
    

    def hasSameUniqueSignatureAs(self, other_isoform):
        return (self.unique_signature == other_isoform.getUniqueSignature())
    

    def setTSSAndSpliceVariantAndPolyASiteNumbers(self, tss_number, splice_variant_number, polyA_site_number):
        self.tss_number = tss_number
        self.splice_variant_number = splice_variant_number
        self.polyA_variant_number = polyA_site_number
        

    def setLengthAndExonTermini(self):
        self.exon_termini_positions = set()
        self.length = 0
        has_1bp_exon = False
        for e in self.ordered_exons:
            exon_start, exon_stop = e.getStartStop()
            if (exon_start == exon_stop):
                has_1bp_exon = True  # See, for instance, isoform ENST00000619410.4 of SCN7A
            self.exon_termini_positions.add(exon_start)
            self.exon_termini_positions.add(exon_stop)
            self.length += e.getLength()
        try:
            assert ((len(self.exon_termini_positions) >= 2 and len(self.exon_termini_positions) % 2 == 0) or has_1bp_exon), \
              "ERROR: In setLengthAndExonTermini() - malformed isoform, problem with number of exon termini"
        except AssertionError as ae:
            pdb.set_trace()

        self.correspondmRNACoordsToGenomicCoords()


    def getTerminiLengthAdjustments(self):
        if (len(self.ordered_exons) == 1):
            first_exon_length_adjustment_5p, last_exon_length_adjustment_3p = self.ordered_exons[0].getLengthAdjustmentAmounts()
        else:
            if (self.strand == '+'):
                first_exon_length_adjustment_5p, first_exon_length_adjustment_3p = self.ordered_exons[0].getLengthAdjustmentAmounts()
                last_exon_length_adjustment_5p, last_exon_length_adjustment_3p = self.ordered_exons[-1].getLengthAdjustmentAmounts()
            else:
                first_exon_length_adjustment_5p, first_exon_length_adjustment_3p = self.ordered_exons[-1].getLengthAdjustmentAmounts()
                last_exon_length_adjustment_5p, last_exon_length_adjustment_3p = self.ordered_exons[0].getLengthAdjustmentAmounts()

            assert (last_exon_length_adjustment_5p == 0)
            assert (first_exon_length_adjustment_3p == 0)

        return (first_exon_length_adjustment_5p, last_exon_length_adjustment_3p)


    def getChromosome(self):
        return self.chromosome


    def getStrand(self):
        return self.strand


    def getDatabaseAndLocus(self):
        return (self.src_db, self.src_db_locus)


    def isFromSameDBAndLocusAs(self, other_isoform):
        other_src_db, other_src_db_locus = other_isoform.getDatabaseAndLocus()
        return (self.src_db_locus != None and self.src_db_locus == other_src_db_locus and self.src_db == other_src_db)
                

    def getSrcDB(self):
        return self.src_db


    def getSrcDBLocus(self):
        return self.src_db_locus


    def addPolyASite(self, pos):
        self.polyA_sites.add(pos)


    def discardPolyASite(self, pos):
        self.polyA_sites.discard(pos)
            

    def getPolyASites(self):
        return self.polyA_sites


    def getNumberOfExons(self):
        return len(self.ordered_exons)


    def getFirstExon(self):
        return self.ordered_exons[0] if (self.strand == '+') else self.ordered_exons[-1]


    def getLastExon(self):
        return self.ordered_exons[0] if (self.strand == '-') else self.ordered_exons[-1]


    def getExon(self, exon_num):
        """Numbering is in the 5'->3' direction. Can be negative."""
        exon = None
        if (self.strand == '+'):
            exon = self.ordered_exons[exon_num]
        else:
            if (exon_num >= 0):
                exon_num = len(ordered_exons)-1-exon_num
                exon = self.ordered_exons[exon_num]
            else:
                exon_num = abs(exon_num)-1
                exon = self.ordered_exons[exon_num]
        return exon

            
    def getStrandwiseStart(self):
        return self.ordered_exons[0].getStart() if (self.strand == '+') else self.ordered_exons[-1].getStop()
        

    def getStrandwiseStop(self):
        return self.ordered_exons[0].getStart() if (self.strand == '-') else self.ordered_exons[-1].getStop()


    def getExonsAsTuples(self, include_strand=False):
        exon_tuples = []
        if (include_strand):
            for e in self.ordered_exons:
                exon_tuples.append( (self.chromosome, e.getStart(), e.getStop(), e.getStrand()) )
        else:
            for e in self.ordered_exons:
                exon_tuples.append( (self.chromosome, e.getStart(), e.getStop()) )

        return exon_tuples
    

    def getExonTuples(self):
        exon_tuples = set()
        for e in self.ordered_exons:
            exon_tuples.add( (e.getStart(), e.getStop()) )
        return exon_tuples


    def getIntronTuples(self):
        intron_tuples = set()
        for i in xrange(1, len(self.ordered_exons)):
            intron_start = self.ordered_exons[i-1].getStop() + 1
            intron_stop = self.ordered_exons[i].getStart() - 1
            intron_tuples.add( (intron_start, intron_stop) )
        return intron_tuples

        
    def getLengthAllExons(self):
        assert (self.length > 0), "Isoform length should be set by now"

        # Confirm
        total_length = sum(map(methodcaller("getLength"), self.ordered_exons))
        assert (total_length == self.length), "Sum of exon lengths does not equal self.length, in RNAIsoform:getLengthAllExons()"

        return total_length

    
    def getExonTerminiPositions(self):
        return self.exon_termini_positions


    def getExonFivePrimeTerminiPositions(self):
        assert (self.strand in ['-','+'])

        ordered_exon_termini_positions = sorted(self.exon_termini_positions)
        if (self.strand == '-'):
            ordered_exon_termini_positions.reverse()

        return ordered_exon_termini_positions[0::2]


    def getExonThreePrimeTerminiPositions(self):
        assert (self.strand in ['-','+'])

        ordered_exon_termini_positions = sorted(self.exon_termini_positions)
        if (self.strand == '-'):
            ordered_exon_termini_positions.reverse()

        return ordered_exon_termini_positions[1::2]


    def getSpliceJunctionPositions(self):
        SJ_positions = set()
        if (len(self.ordered_exons) > 1):
            ordered_exon_termini_positions = sorted(list(self.exon_termini_positions))
            SJ_positions.update(ordered_exon_termini_positions[1:-1])
        return SJ_positions


    def hasSameSpliceJunctionPositionsAs(self, other_isoform):
        return (self.getSpliceJunctionPositions() == other_isoform.getSpliceJunctionPositions())

        
    def getStartStop(self):
        ordered_exon_termini_positions = sorted(list(self.exon_termini_positions))
        return (ordered_exon_termini_positions[0], ordered_exon_termini_positions[-1])

        
    def getTranscriptIDInSrcDB(self):
        return self.ID

    
    def getLocusNames(self):
        return self.locus_names


    def getLocusSymbolFromDatabase(self, src_db):
        assert (src_db in ["RefSeq", "GENCODE", "HInv", "AceView", "SIB", "ALTSCAN", "UCSC"])
        locus_symbol = None
        
        if (src_db in self.locus_names):
            locus_symbol = self.locus_names[src_db][0]
              
        return locus_symbol


    def getLocusSymbol(self):
        """Return the preferred non-CGDB name of this isoform """
        locus_symbol = None

        if (self.hgnc_symbol != None):
            locus_symbol = self.hgnc_symbol
        elif ("RefSeq" in self.locus_names):
            locus_symbol = self.locus_names["RefSeq"][0]
        elif ("GENCODE" in self.locus_names):
            locus_symbol = self.locus_names["GENCODE"][0]
        elif ("HInv" in self.locus_names):
            locus_symbol = self.locus_names["HInv"][0]
        elif ("AceView" in self.locus_names):
            locus_symbol = self.locus_names["AceView"][0]
              
        return locus_symbol

        
    def getHGNCID(self):
        return self.hgnc_id


    def getHGNCSymbol(self):
        return self.hgnc_symbol


    def setHGNCID(self, hgnc_id):
        assert (self.hgnc_id == None)
        self.hgnc_id = hgnc_id


    def setHGNC(self, hgnc_id, hgnc_symbol):
        assert (self.hgnc_id == None and self.hgnc_symbol == None)
        self.hgnc_id = hgnc_id
        self.hgnc_symbol = hgnc_symbol


    def updateHGNC(self, hgnc_id, hgnc_symbol):
        if (self.hgnc_id == None):
            self.hgnc_id = hgnc_id
        else:
            assert (self.hgnc_id == hgnc_id)
        
        if (self.hgnc_symbol == None):
            self.hgnc_symbol = hgnc_symbol
        else:
            assert (self.hgnc_symbol == hgnc_symbol)

        
    def setCGDBName(self, isoform_group_number):
        assert (self.tss_number != None), "Setting CGDB name, missing TSS number"
        assert (self.splice_variant_number != None), "Setting CGDB name, missing splice variant number"
        assert (self.polyA_variant_number != None), "Setting CGDB name, missing polyA variant number"
        assert (self.CGDB_name == None), "Setting CGDB name, already set"
        self.CGDB_name = "CG_%d.%d.%d.%d" % \
          (isoform_group_number, self.tss_number, self.splice_variant_number, self.polyA_variant_number)
        
        
    def hasCGDBName(self):
        return self.CGDB_name != None


    def getCGDBName(self):
        name = self.CGDB_name if (self.CGDB_name != None) else self.ID
        assert (name.startswith("CG_"))
        return name

    
    def setFivePrimeIncomplete(self):
        self.is_five_prime_complete = False
        

    def setThreePrimeIncomplete(self):
        self.is_three_prime_complete = False


    def setFivePrimeComplete(self):
        self.is_five_prime_complete = True
        

    def setThreePrimeComplete(self):
        self.is_three_prime_complete = True


    def isFivePrimeComplete(self):
        return self.is_five_prime_complete


    def isThreePrimeComplete(self):
        return self.is_three_prime_complete


    def addLocusNames(self, database, names):
        if (database not in self.locus_names):
            self.locus_names[database] = []

        if (isinstance(names,list)):
            self.locus_names[database].extend(alt_names_list)
        else:
            self.locus_names[database].append(names)

            
    def addExons(self, chromosome, start, block_sizes, block_starts, strand):
        assert (chromosome == self.chromosome)
        assert (strand == self.strand)
        assert (isinstance(block_sizes, list) and isinstance(block_starts, list))
            
        for rel_exon_start, exon_len in zip(block_starts, block_sizes):
            abs_exon_start = start + rel_exon_start 
            abs_exon_stop = abs_exon_start + exon_len
            self.addExon(chromosome, abs_exon_start+1, abs_exon_stop, strand) # "+1" because BED start positions are 0-based


    def addExon(self, chromosome, start, stop, strand):
        """Maintains exons in their strandwise physical order. For a '-'/'+' strand
        isoform, the first exon in self.ordered_exons will be the last/first exon.
        Instantiated after all exons read by call to getFirstExon() in setIsoformCompleteness()"""
        assert (self.chromosome == chromosome), "Inconsistent chromosome"
        assert (self.strand == strand), "Inconsistent exon strand"
        
        new_exon = Exon(self, chromosome, start, stop, strand)

        self.ordered_exons.append(new_exon)
        self.ordered_exons = sorted(self.ordered_exons, key=attrgetter("start"))


    def getIntronicSignature(self):
        """Returns a tuple of the first and last intron positions in sorted order, regardless
        of the strand of the mRNA (i.e. - or + strand). The strand is added to the tuple if requested.
        """
        intronic_signature = None
        if (len(self.ordered_exons) > 1):
            positions = [self.ordered_exons[0].getStop(), self.ordered_exons[-1].getStart()]
            for exon in self.ordered_exons[1:-1]:
                positions.extend( [exon.getStart(), exon.getStop()] )
            positions.sort()
            intronic_signature = tuple(positions)

        return intronic_signature


    def getStrandwiseSequence(self, num_nucs_needed, genome_seq): #, as_string=False):
        """num_nucs_needed is either 'all' or an integer"""
        nuc_seq = "" #Seq("", IUPAC.unambiguous_dna)

        exons_5p_3p = self.ordered_exons
        if (self.strand == "-"):
            exons_5p_3p = exons_5p_3p[::-1]
            
        for exon in exons_5p_3p:
            nuc_seq += exon.getStrandwiseSequence(genome_seq) 
            if (num_nucs_needed != 'all' and len(nuc_seq) >= num_nucs_needed):
                break
    
        #if (as_string):
        #    nuc_seq = str(nuc_seq)

        return nuc_seq


    def setSequenceFromGenome(self, genome_seq):
        self.nuc_seq = self.getStrandwiseSequence("all", genome_seq) #, True)
    

    def setSequenceFromSequence(self, nuc_seq):
        assert (self.length == len(nuc_seq)), "Nucleotide sequence not equal to sum of all exon lengths"
        self.nuc_seq = nuc_seq


    def getThreePrimeTerminalRegionSequence(self, genome_seq, flanking_region_size):
        """Returns the sequence region centered on the last nucleotide of the 3' exon. Sequence returned is in 
        the 5'->3' orientation. That is, if isoform is on the negative strand, the sequence will be reverse complemented. """
        terminal_region_nuc_seq = None
        terminus_three_prime = self.getStrandwiseStop()

        needed_len = flanking_region_size + 1
        length_all_exons = self.getLengthAllExons()
        if (length_all_exons >= needed_len):
            terminal_region_nuc_seq = ""
            accum_len = 0
            i = -1
            while (accum_len < needed_len):
                exon = self.getExon(i)
                accum_len += exon.getLength()
                exon_seq = exon.getStrandwiseSequence(genome_seq, True)
                terminal_region_nuc_seq = exon_seq + terminal_region_nuc_seq
                assert (accum_len == len(terminal_region_nuc_seq))
            terminal_region_nuc_seq = terminal_region_nuc_seq[-needed_len:]    

            #if (self.strand == '+'):
            #    genome_region_spec = "%s:%d-%d" % (self.chromosome, terminus_three_prime+1, terminus_three_prime+flanking_region_size)
            #else:
            #    genome_region_spec = "%s:%d-%d" % (self.chromosome, terminus_three_prime-flanking_region_size, terminus_three_prime-1)
            #nuc_seq_fasta = pysam.faidx(genome_fasta, genome_region_spec)
            #nuc_seq = ''.join(map(lambda x: x.strip(), nuc_seq_fasta[1:]))

            if (self.strand == '+'):
                nuc_seq = genome_seq[self.chromosome][terminus_three_prime:terminus_three_prime+flanking_region_size].seq.upper()
            else:
                nuc_seq = genome_seq[self.chromosome][terminus_three_prime-flanking_region_size-1:terminus_three_prime-1].seq.upper()

            if (self.strand == '-'):
                nuc_seq = Seq(nuc_seq, IUPAC.unambiguous_dna)
                nuc_seq = nuc_seq.reverse_complement()
                nuc_seq = nuc_seq.tostring()
                
            terminal_region_nuc_seq += nuc_seq.upper()

        return terminal_region_nuc_seq


    def getFirstExonForLabeledBED(self):
        first_exon = self.getFirstExon()
        exon_start = first_exon.getStart() - 1
        exon_stop = first_exon.getStop()
        this_ID = "exon|%s|%d-%d|%s" % (self.ID, exon_start, exon_stop, self.strand)
        bed_line = "%s\t%d\t%d\t%s\t0\t%s" % (self.chromosome, exon_start, exon_stop, this_ID, self.strand)
        return (this_ID, bed_line)
        

    def getLastExonForLabeledBED(self):
        last_exon = self.getLastExon()
        exon_start = last_exon.getStart() - 1
        exon_stop = last_exon.getStop()
        this_ID = "exon|%s|%d-%d|%s" % (self.ID, exon_start, exon_stop, self.strand)
        bed_line = "%s\t%d\t%d\t%s\t0\t%s" % (self.chromosome, exon_start, exon_stop, this_ID, self.strand)
        return (this_ID, bed_line)
        

    def getFirstExonStrandwiseStartForLabeledBED(self):
        first_exon = self.getFirstExon()
        exon_start = first_exon.getStart() - 1
        exon_stop = first_exon.getStop()
        this_ID = "exon|%s|%d-%d|%s" % (self.ID, exon_start, exon_stop, self.strand)

        stop = first_exon.getStrandwiseStart()
        start = stop - 1
        bed_line = "%s\t%d\t%d\t%s\t0\t%s" % (self.chromosome, start, stop, this_ID, self.strand)

        return (this_ID, bed_line)


    def getLastExonStrandwiseStopForLabeledBED(self):
        last_exon = self.getLastExon()
        exon_start = last_exon.getStart() - 1
        exon_stop = last_exon.getStop()
        this_ID = "exon|%s|%d-%d|%s" % (self.ID, exon_start, exon_stop, self.strand)
        
        stop = last_exon.getStrandwiseStop()
        start = stop - 1
        bed_line = "%s\t%d\t%d\t%s\t0\t%s" % (self.chromosome, start, stop, this_ID, self.strand)

        return (this_ID, bed_line)


    def getAsBED12Fields(self, naming_scheme, CDS=None, color_triple="0,0,0"):
        all_exon_tups = map(lambda x: x.getStartStop(), self.ordered_exons)
        all_exon_tups = sorted(all_exon_tups, key=lambda x:x[0])

        chromStart = all_exon_tups[0][0]
        chromEnd = all_exon_tups[-1][1]
        blockCount = len(all_exon_tups)

        if (CDS == None):
            CDS = (chromStart-1, chromStart-1)

        blockStarts = []
        blockSizes = []
        for s,e in all_exon_tups:
            blockStarts.append( s - chromStart )
            blockSizes.append( e - s + 1 )

        if (naming_scheme == "CGDB"):
            isoform_name = self.CGDB_name
        elif (naming_scheme == "ID"):
            isoform_name = self.ID
        elif (naming_scheme == "locus_and_ID"):
            isoform_name = "%s:%s" % (self.src_db, self.ID)
        elif (naming_scheme == "unique_signature"):
            isoform_name = "%s:%s:%s" % (self.src_db, self.ID, self.unique_signature)
        elif (naming_scheme == "instance_address"):
            isoform_name = str(self)
        else:
            print("ERROR: unrecognized isoform naming scheme \'%s\' for BED12. Exiting" % naming_scheme, file=sys.stderr)
            sys.exit(1)
            
        bed_fields = [self.chromosome, str(chromStart-1), str(chromEnd), isoform_name, "0", self.strand,
                     str(CDS[0]), str(CDS[1]), color_triple, str(blockCount),
                     ",".join(map(str, blockSizes)), ",".join(map(str, blockStarts))] #  For debugging , self.src_db, str(self.ID)

        return bed_fields


    def getAsGTF(self, naming_scheme, score=1000):
        lines = []

        if (naming_scheme == "CGDB"):
            isoform_name = self.CGDB_name
        elif (naming_scheme == "ID"):
            isoform_name = self.ID
        elif (naming_scheme == "locus_and_ID"):
            isoform_name = "%s:%s" % (self.src_db, self.ID)
        elif (naming_scheme == "unique_signature"):
            isoform_name = "%s:%s:%s" % (self.src_db, self.ID, self.unique_signature)
        else:
            print("ERROR: unrecognized isoform naming scheme \'%s\' for BED12. Exiting" % naming_scheme, file=sys.stderr)
            sys.exit(1)

        #
        # First is the bare and unextended exon structure
        #
        for exon in self.ordered_exons:
            if (exon.isFirstExon()):
                if (self.strand == '-'):
                    s = exon.getStart()
                    e = exon.getOriginalStrandwiseStart()
                else:
                    s = exon.getOriginalStrandwiseStart()
                    e = exon.getStop()
            else:
                s,e = exon.getStartStop()

            gtf_elems = [self.chromosome, "CGDB", "exon", str(s), str(e), str(score), self.strand,
                        ".", "gene_id \"%s\"; transcript_id \"%s\";" % (self.src_db_locus, isoform_name)]
            lines.append( "\t".join(gtf_elems) )

        #
        # Next is any special annotation to add on top of the bare exon structure
        #
        
        # The extended 1st exon, if applicable
        first_exon = self.getFirstExon()
        curr_start = first_exon.getStrandwiseStart()
        orig_start = first_exon.getOriginalStrandwiseStart()
        if (curr_start != orig_start):
            if (self.strand == '-'):
                s = min(curr_start, orig_start) + 1
                e = max(curr_start, orig_start)
            else:
                s = min(curr_start, orig_start)
                e = max(curr_start, orig_start) - 1
            gtf_elems = [self.chromosome, "CGDB", "CDS", str(s), str(e), str(score), self.strand, ".",
                         "gene_id \"%s\"; transcript_id \"%s\";" % (self.src_db_locus, isoform_name)]
            lines.append( "\t".join(gtf_elems) )

        return lines


    def getAsGFF3(self, parent_gff3_ID):
        mRNA_tuple, exon_tuples = None, []

        start, stop = self.getStartStop()
        CGDB_name = self.getCGDBName()
        assert (CGDB_name != None)
        if (parent_gff3_ID != None):
            mRNA_tuple = (self.chromosome, start, stop, self.strand, "ID=%s;Parent=%s;" % (CGDB_name, parent_gff3_ID))
        else:
            mRNA_tuple = (self.chromosome, start, stop, self.strand, "ID=%s;" % CGDB_name)

        for exon in self.ordered_exons:
            exon_tuple = exon.getAsGFF3()
            exon_tuples.append( exon_tuple )
                        
        return (mRNA_tuple, exon_tuples)


    def correspondmRNACoordsToGenomicCoords(self, set_rna_to_genome=True, set_genome_to_rna=False):
        '''Creates a 0-based mRNA index to 1-based genomic index tuple list.'''
        mRNA_indices = range(self.length)

        genomic_indices = []
        for exon in self.ordered_exons:
            exon_start, exon_stop = exon.getStartStop()
            genomic_indices.extend( range(exon_start, exon_stop+1) )

        assert (self.length != 0 and self.length != None)
        assert (len(genomic_indices) == self.length)
            
        self.genomic_positions = set(genomic_indices)
        
        genomic_indices = sorted(genomic_indices, reverse=(self.strand=='-'))
        if (set_rna_to_genome):
            self.coords_mRNA_to_genome = np.array(genomic_indices, dtype=np.uint32)
        if (set_genome_to_rna):
            self.coords_genome_to_mRNA = dict( zip(genomic_indices, range(self.length)) )


    def clearCorrespondmRNACoordsToGenomicCoords(self, clear_rna_to_genome=False, clear_genome_to_rna=True):
        self.genomic_positions = None
        if (clear_rna_to_genome):
            self.coords_mRNA_to_genome = None
        if (clear_genome_to_rna):
            self.coords_genome_to_mRNA = None


    def getFirstExonGenomicPositions(self):
        exon = self.ordered_exons[0] if (self.strand == '+') else self.ordered_exons[-1]
        start, stop = exon.getStartStop()
        genomic_positions = set(range(start,stop+1))
        return genomic_positions


    def getLastExonGenomicPositions(self):
        genomic_positions = set()
        exon = self.ordered_exons[-1] if (self.strand == '+') else self.ordered_exons[0]
        start, stop = exon.getStartStop()
        genomic_positions = set(range(start,stop+1))
        return genomic_positions


    def getExonsGenomicPositions(self):
        genomic_positions_per_exon = []

        pos = 0
        if (self.strand == '+'):
            for e in self.ordered_exons:
                exon_len = e.getLength()
                genomic_positions_per_exon.append( set(map(lambda i: self.coords_mRNA_to_genome[i], range(pos, pos+exon_len))) )
                pos += exon_len
        else:
            for e in self.ordered_exons[::-1]:
                exon_len = e.getLength()
                genomic_positions_per_exon.append( set(map(lambda i: self.coords_mRNA_to_genome[i], range(pos, pos+exon_len))) )
                pos += exon_len

        return genomic_positions_per_exon


    def getGenomicPositions5p3p(self):
        return self.coords_mRNA_to_genome


    def getGenomicCoord(self, mRNA_pos):
        return self.coords_mRNA_to_genome[mRNA_pos]


    def convertToGenomicCoords(self, mRNA_positions):
        return list(map(lambda mRNA_pos: self.coords_mRNA_to_genome[mRNA_pos], mRNA_positions))


    def convertToTranscriptCoords(self, genomic_positions):
        return list(map(lambda pos: self.coords_genome_to_mRNA[pos], genomic_positions))


    def getAllGenomicCoords(self):
        return set(self.coords_genome_to_mRNA.keys())


    def getIsoformCoord(self, genomic_pos):
        return self.coords_genome_to_mRNA[genomic_pos]


    def getmRNAPositions(self, genomic_position_5p, genomic_position_3p):
        mRNA_position_5p = self.coords_genome_to_mRNA[genomic_position_5p]
        mRNA_position_3p = self.coords_genome_to_mRNA[genomic_position_3p]
        return (mRNA_position_5p, mRNA_position_3p)


    def getSequenceLength(self, genomic_position_5p, genomic_position_3p):
        if (genomic_position_5p in self.coords_genome_to_mRNA and genomic_position_3p in self.coords_genome_to_mRNA):
            mRNA_position_5p = self.coords_genome_to_mRNA[genomic_position_5p]
            mRNA_position_3p = self.coords_genome_to_mRNA[genomic_position_3p]
            assert (mRNA_position_5p <= mRNA_position_3p), "mRNA 5' position > 3' position"
            seqlen = mRNA_position_3p - mRNA_position_5p + 1
        else:
            seqlen = None
        return seqlen


    def getNegativeSequenceLength(self, genomic_position_5p, genomic_position_3p):
        '''Returns the length of sequence on the mRNA between including the transcript positions corresponding to the genomic positions'''
        if (genomic_position_5p in self.coords_genome_to_mRNA and genomic_position_3p in self.coords_genome_to_mRNA):
            mRNA_position_5p = self.coords_genome_to_mRNA[genomic_position_5p]
            mRNA_position_3p = self.coords_genome_to_mRNA[genomic_position_3p]
            assert (mRNA_position_5p > mRNA_position_3p)
            neg_seqlen = mRNA_position_3p - mRNA_position_5p - 1
        else:
            neg_seqlen = None
        return neg_seqlen


    def getSequence(self, genome_seq, genomic_position_5p, genomic_position_3p):
        mRNA_position_5p = self.coords_genome_to_mRNA[genomic_position_5p]
        mRNA_position_3p = self.coords_genome_to_mRNA[genomic_position_3p]

        assert (mRNA_position_5p <= mRNA_position_3p), "mRNA 5' position > 3' position"
        if (self.nuc_seq == None):
            self.setSequence(genome_seq)

        return self.nuc_seq[mRNA_position_5p:mRNA_position_3p + 1]


    def getValidSubsequence(self, genomic_position_1, genomic_position_2):
        assert (genomic_position_1 != genomic_position_2)
        mRNA_position_5p = self.coords_genome_to_mRNA[genomic_position_1]
        mRNA_position_3p = self.coords_genome_to_mRNA[genomic_position_2]
        if (mRNA_position_5p < mRNA_position_3p):
            subseq = self.nuc_seq[mRNA_position_1:mRNA_position_2 + 1]
        else:
            subseq = self.nuc_seq[mRNA_position_2:mRNA_position_1 + 1]
        return subseq


    # This method avoids the complexities of getPCRProduct() with multiple primer matches in the template
    def hasPCRProduct(self, fwd_primer_seq, rev_primer_seq, amplicon, allow_primer_mismatches=False):
        try:
            rev_primer_seq_revcomp = translate(rev_primer_seq, DNA_complement_table)[::-1]
        except NameError:  # TODO: Remove when all upgraded to python3
            import string
            rev_primer_seq_revcomp = rev_primer_seq[::-1].translate(DNA_complement_table)

        exact_amplicon_found = amplicon in self.nuc_seq and amplicon.startswith(fwd_primer_seq) and amplicon.endswith(rev_primer_seq_revcomp)
        if (not exact_amplicon_found and allow_primer_mismatches):
            inexact_amplicon_found = amplicon[len(fwd_primer_seq):-len(rev_primer_seq)] in self.nuc_seq
        else:
            inexact_amplicon_found = False

        return exact_amplicon_found or inexact_amplicon_found


    def getPCRProductTranscriptCoords(self, fwd_primer_seq, rev_primer_seq, amplicon, allow_primer_mismatches=False):
        try:
            assert( self.hasPCRProduct(fwd_primer_seq, rev_primer_seq, amplicon, allow_primer_mismatches) )
        except AssertionError:
            pdb.set_trace()
            self.hasPCRProduct(fwd_primer_seq, rev_primer_seq, amplicon, allow_primer_mismatches)
        insert = amplicon[len(fwd_primer_seq):-len(rev_primer_seq)]
        insert_start_index = self.nuc_seq.index(insert)
        return list(range(insert_start_index-len(fwd_primer_seq),insert_start_index+len(insert)+len(rev_primer_seq)))


    def getPCRProduct(self, fwd_primer_seq, rev_primer_seq):
        '''The rev_primer_seq is the actual reverse primer sequence, which is the reverse complement of the mRNA subsequence.'''
        print("TODO: A check needs to be put in place against having multiple fwd and/or rev primer matches. Maybe base on thermodynamic frac_duplexed.", file=sys.stderr)
        sys.exit(1)
        primers_product = None
        rev_primer_seq_revcomp = translate(rev_primer_seq, DNA_complement_table)[::-1]
        fwd_primer_start = self.nuc_seq.find(fwd_primer_seq)
        rev_primer_start = self.nuc_seq.find(rev_primer_seq_revcomp)
        if (fwd_primer_start != -1 and rev_primer_start != -1):
            rev_primer_end = rev_primer_start + len(rev_primer_seq_revcomp)
            primers_product = self.nuc_seq[fwd_primer_start:rev_primer_end]
        return primers_product


    # TODO: A check needs to be put in place against having multiple fwd and/or rev primer matches.
#    def getPCRProductBED12(self, fwd_primer_seq, rev_primer_seq, chromosome, label, score=0, strand = ".", itemRgb="0,0,0"):
#        '''The rev_primer_seq is the actual reverse primer sequence, which is the reverse complement of the mRNA subsequence.'''
#        bed12_line = None
#        rev_primer_seq_revcomp = translate(rev_primer_seq, DNA_complement_table)[::-1]
#        fwd_primer_start = self.nuc_seq.find(fwd_primer_seq)
#        rev_primer_start = self.nuc_seq.find(rev_primer_seq_revcomp)
#        if (fwd_primer_start != -1 and rev_primer_start != -1):
#            rev_primer_end = rev_primer_start + len(rev_primer_seq_revcomp)
#            primers_product_genomic_positions = self.coords_mRNA_to_genome[fwd_primer_start:rev_primer_end]
#            primers_product_genomic_positions.sort()
#
#            start = primers_product_genomic_positions[0]
#            stop = primers_product_genomic_positions[-1]
#
#            block_sizes = []
#            block_starts = []
#            for k, g in groupby(enumerate(primers_product_genomic_positions), lambda (i,x):i-x):
#                    group = map(itemgetter(1), g)
#                    block_starts.append( str(group[0] - start) )
#                    block_sizes.append( str(group[-1] - group[0] + 1) )
#        
#            bed12_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s" % \
#                                 (chromosome, start-1, stop, label, score, strand, start, start, itemRgb, len(block_starts), ','.join(block_sizes), ','.join(block_starts))
#        return bed12_line


    def numSharedPositions(self, genomic_positions):
        return len(self.genomic_positions.intersection(genomic_positions))


    def hasPositions(self, genomic_positions):
        return genomic_positions.issubset(self.genomic_positions)
