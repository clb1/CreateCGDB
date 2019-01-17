class TSS:
    def __init__(self, source_db, source_label, tss_label, chromosome, start, stop, strand):
        """Expects start and stop to be in 1-based coordinates."""
        self.src_db = source_db
        self.src_label = source_label
        self.ID = tss_label
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.halfwidth = (stop - start)/2
        self.strand = strand

    def isSrcLabel(self, src_db, src_label):
        return self.src_db == src_db and self.src_label == src_label
    
    def getSrcDB(self):
        return self.src_db

    def getSrcLabel(self):
        return self.src_label
    
    def getChromosome(self):
        return self.chromosome

    def getStrand(self):
        return self.strand

    def getStart(self):
        return self.start

    def getStop(self):
        return self.stop

    def getFivePrimePosition(self):
        return self.start if (self.strand == '+') else self.stop

    def getThreePrimePosition(self):
        return self.stop if (self.strand == '-') else self.start

    def getFivePrimeThreePrimePositions(self):
        return (self.start, self.stop) if (self.strand == '+') else (self.stop, self.start)
    
    def doesOverlap(self, pos):
        return pos >= self.start and pos <= self.stop
         
    def isWithinHalfwidth(self, pos):
        """Returns a triple of 1) absolute distance of pos from an edge of the TSS,
        2) the relative direction of the TSS from pos, and 3) a reference to this TSS."""
        ret_val = None
        if (pos >= self.start and pos <= self.stop):
            ret_val = (0, "NA", self)
        elif (pos < self.start and self.start - pos <= self.halfwidth):
            ret_val = (self.start - pos, "downstream" if (self.strand == '+') else "upstream", self)
        elif (pos > self.stop and pos - self.stop <= self.halfwidth):
            ret_val = (pos - self.stop, "upstream" if (self.strand == '+') else "downstream", self)
        return ret_val

    def getTSSForLabeledBED(self, slop=0.0):
        this_ID = "TSS|%s|%s|%s" % (self.src_db, self.src_label, self.ID)
        if (slop==0.0):
            bed_line = "%s\t%d\t%d\t%s\t0\t%s" % (self.chromosome, self.start-1, self.stop, this_ID, self.strand)
        else:
            ext_len = int(slop * (self.stop - self.start + 1))
            bed_line = "%s\t%d\t%d\t%s\t0\t%s" % (self.chromosome, self.start-1-ext_len, self.stop+ext_len, this_ID, self.strand)
        return (this_ID, bed_line)
                
