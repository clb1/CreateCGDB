import sys

class PolyA(object):
    def __init__(self, source_db, pAs_label, chromosome, start, stop, strand, cluster_median, cluster_mode):
        """Expects start and stop to be in 1-based coordinates."""
        self.src_db = source_db
        self.ID = pAs_label
        self.chromosome = chromosome
        self.strand = strand
        self.start = start
        self.stop = stop
        self.cluster_median = cluster_median
        self.cluster_mode = cluster_mode

        #self.halfwidth = (stop - start)/2


    def isSrcDB(self, src_db):
        return self.src_db == src_db
    
    def getSrcDB(self):
        return self.src_db

    def getChromosome(self):
        return self.chromosome

    def getStrand(self):
        return self.strand

    def getStart(self):
        return self.start

    def getStop(self):
        return self.stop

    def getSitePosition(self):
        if (self.src_db == "APADB"):
            return self.cluster_mode
        else:
            print >> sys.stderr, "ERROR: have not implemented polyA site definition for database %s. Exiting." % self.src_db
            sys.exit(1)
            
    def getFivePrimePosition(self):
        return self.start if (self.strand == '+') else self.stop

    def getThreePrimePosition(self):
        return self.stop if (self.strand == '-') else self.start

    def getFivePrimeThreePrimePositions(self):
        return (self.start, self.stop) if (self.strand == '+') else (self.stop, self.start)
    
    def doesOverlap(self, pos):
        return pos >= self.start and pos <= self.stop
         
#    def isWithinHalfwidth(self, pos):
#        """Returns a triple of 1) absolute distance of pos from an edge of the TSS,
#        2) the relative direction of the TSS from pos, and 3) a reference to this TSS."""
#        print >> sys.sterr, "ERROR: halfwidth usage not yet resolved."
#        sys.exit(1)
#        ret_val = None
#        if (pos >= self.start and pos <= self.stop):
#            ret_val = (0, "NA", self)
#        elif (pos < self.start and self.start - pos <= self.halfwidth):
#            ret_val = (self.start - pos, "downstream" if (self.strand == '+') else "upstream", self)
#        elif (pos > self.stop and pos - self.stop <= self.halfwidth):
#            ret_val = (pos - self.stop, "upstream" if (self.strand == '+') else "downstream", self)
#        return ret_val

    def getPolyASiteForLabeledBED(self, slop=0):
        this_ID = "pAs|%s|%s" % (self.src_db, self.ID)

        #if (slop==0.0):
        #    bed_line = "%s\t%d\t%d\t%s\t%d\t%s" % (self.chromosome, self.cluster-1, self.stop, this_ID, self.cluster_mode, self.strand)
        #else:
        bed_line = "%s\t%d\t%d\t%s\t%d\t%s" % (self.chromosome, self.cluster_mode-1-slop, self.cluster_mode+slop, this_ID, 0, self.strand)
        return (this_ID, bed_line)
                
