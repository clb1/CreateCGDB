from collections import defaultdict
from itertools import chain
import gzip
import bz2
import os
import re
import sys
from RNAIsoform_python2 import RNAIsoform
from RNAIsoformAnnotator import RNAIsoformAnnotator
from TSS import TSS
from PolyA import PolyA

import pdb

def readIsoformModels(isoform_model_databases, annotator):
    all_isoform_models = defaultdict(set)
    
    for db in isoform_model_databases.split(','):
        print >> sys.stderr, "\nINFO: parsing %s" % db
        isoform_models = parseDatabase(db, annotator)
        for isoform in isoform_models:
            okay_to_use = isoform.setIsoformCompleteness()
            if (okay_to_use):
                chromosome = isoform.getChromosome()
                all_isoform_models[chromosome].add( isoform )
            else:
                print >> sys.stderr, "WARNING: problem with completeness of %s. Skipping." % isoform.getTranscriptIDInSrcDB()

    print >> sys.stderr, "Checking for strand inconsistencies between isoforms from different databases"
    checkAndCorrectStrandInconsistencies(all_isoform_models)

    return all_isoform_models


def checkAndCorrectStrandInconsistencies(all_isoform_models):
    isoforms_to_remove = defaultdict(list)
    for chromosome, chrom_isoforms in all_isoform_models.items():
        if (len(chromosome) > 5):
            continue
        print >> sys.stderr, "# %s #" % chromosome
        isoform_w_same_intronic_sig = defaultdict(list)
        for isoform in chrom_isoforms:
            sig = isoform.getIntronicSignature()
            if (sig != None): # Skip single-exon isoforms
                isoform_w_same_intronic_sig[sig].append(isoform)

        for sig, sig_isoforms in isoform_w_same_intronic_sig.items():
            if (len(set(map(lambda x: x.getStrand(), sig_isoforms))) == 2): # Inconsistent strands for same intronic signature
                #print >> sys.stderr, "%s\t%s" % (str(sig), " ".join(map(lambda x: x.getTranscriptIDInSrcDB(), sig_isoforms)))
                if (len(sig_isoforms)==2):
                    DBs = set(map(lambda i: i.getSrcDB(), sig_isoforms))
                    if (len(DBs) == 1):
                        DB = list(DBs)[0]
                        isoforms_to_remove[chromosome].extend(sig_isoforms) # Inconsistency internal to one DB. No easy way to resolve, so just eliminate problem.
                        print >> sys.stderr, "Removing %s isoforms %s" % (DB, " ".join(map(lambda x: x.getTranscriptIDInSrcDB(), sig_isoforms)))
                    else:
                        D = dict(map(lambda i: (i.getSrcDB(),i.getStrand()), sig_isoforms))
                        for src_db in ["RefSeq", "GENCODE", "AceView", "SIB", "ALTSCAN", "RNAcentral", "UCSC", "HInv"]: # Ranked by (empirically observed) trustworthiness
                            if (D.has_key(src_db)):
                                correct_strand = D[src_db]
                                break
                        for sig_isoform in filter(lambda i: i.getStrand() != correct_strand, sig_isoforms):
                            print >> sys.stderr, "Switching strand for %s" % sig_isoform.getTranscriptIDInSrcDB()
                            sig_isoform.switchStrand()
                else:
                    strands = map(lambda i: i.getStrand(), sig_isoforms)
                    num_plus = strands.count('+')
                    num_minus = strands.count('-')
                    if (num_plus != num_minus): # If there is a majority, go with that. Otherwise, leave as is.
                        majority_strand = '+' if (num_plus > num_minus) else '-'
                        for sig_isoform in filter(lambda i: i.getStrand() != majority_strand, sig_isoforms):
                            print >> sys.stderr, "Switching strand for %s" % sig_isoform.getTranscriptIDInSrcDB()
                            sig_isoform.switchStrand()

    for chromosome, isoforms in isoforms_to_remove.items():
        for isoform in isoforms:
            all_isoform_models[chromosome].remove(isoform)



def readTSSs(tss_databases):
    all_TSS = defaultdict(list)

    for db in tss_databases.split(','):
        print >> sys.stderr, "\nINFO: parsing %s" % db
        tss = parseDatabase(db)
        for a_tss in tss:
            chromosome = a_tss.getChromosome()
            all_TSS[chromosome].append(a_tss)

    return all_TSS


def readPolyASites(polyAsite_databases):
    all_polyAsites = defaultdict(list)

    for db in polyAsite_databases.split(','):
        print >> sys.stderr, "\nINFO: parsing %s" % db
        pAs = parseDatabase(db)
        for a_pAs in pAs:
            chromosome = a_pAs.getChromosome()
            all_polyAsites[chromosome].append(a_pAs)

    return all_polyAsites


def parseDatabase(db, annotator=None):
    db_models = None
    db_name = None
    
    if (os.path.basename(db).startswith("AceView")):
        db_name = "AceView"
        db_models = parseAceView(db, annotator)
    elif (os.path.basename(db).startswith("GENCODE")):
        db_name = "GENCODE"
        db_models = parseGENCODE(db, annotator) 
    elif (os.path.basename(db).startswith("UCSC")):
        db_name = "UCSC"
        db_models = parseUCSC(db, annotator)
    elif (os.path.basename(db).startswith("ALTSCAN")):
        db_name = "ALTSCAN"
        db_models = parseALTSCAN(db, annotator)
    elif (os.path.basename(db).startswith("RefSeq")):
        db_name = "RefSeq"
        db_models = parseRefSeq(db, annotator)
    elif (os.path.basename(db).startswith("SIB")):
        db_name = "SIB"
        db_models = parseSIB(db, annotator)
    elif (os.path.basename(db).startswith("HInv")):
        db_name = "HInv"
        db_models = parseHInv(db, annotator)
    elif (os.path.basename(db).startswith("RNAcentral")):
        db_name = "RNAcentral"
        db_models = parseRNAcentral(db)
    elif (os.path.basename(db).startswith("CGDB")):
        db_name = "CGDB"
        db_models = parseCGDB(db, annotator)
    elif (db.startswith("FANTOM5")):
        db_name = "FANTOM5"
        db_models = parseFANTOM(db)
    elif (db.startswith("APADB")):
        db_name = "APADB"
        db_models = parseAPADB(db)
    else:
        print >> sys.stderr, "ERROR: no parser for database %s. Exiting." % db
        sys.exit(1)
        
    print >> sys.stderr, "INFO: read %d models from %s database" % (len(db_models), db_name)

    return db_models


def openDatabaseFile(db):
    ip = None
    if (db.endswith(".gz")):
        ip = gzip.open(db, "rb")
    elif (db.endswith(".bz2")):
        ip = bz2.BZ2File(db, 'r')
    else:
        ip = open(db, 'r')
    return ip


def parseAceView(db, annotator):
    isoform_models = {}
    ip = openDatabaseFile(db)

    loci_to_remove = set()
    
    reExonAnnotation = re.compile("gene_id (\S+).\s+.+\s+transcript_id (\S+). exon_number \d+")
    reTranscriptAnnotation = re.compile("gene_id (\S+). transcript_id (\S+). mRNA_start_support (\S). mRNA_end_support (\S). number_exons \d+")
    reGENCODE = re.compile("GENCODE (\S+)")
    reRefSeq = re.compile("RefSeq (\S+)")

    for line in ip:
        fields = line.strip().split("\t")
        if (fields[2] == "mRNA"):
            chromosome = fields[0]
            strand = fields[6]
            mo = reTranscriptAnnotation.match(fields[-1])
            assert (mo != None)
            gene_id, transcript_id, start_supported, end_supported = mo.groups()

            assert (not isoform_models.has_key(transcript_id)), "Aceview transcript %s already created" % transcript_id
            new_isoform = RNAIsoform("AceView", gene_id, transcript_id, chromosome, strand)

            added_gencode_or_refseq = False
            mo = reGENCODE.search(fields[-1])
            if (mo != None):
                gencode_id = mo.group(1)
                gencode_id = gencode_id.replace(';','')
                new_isoform.addLocusNames("GENCODE", gencode_id)
                added_gencode_or_refseq = True

            mo = reRefSeq.search(fields[-1])
            if (mo != None):
                refseq_id = mo.group(1)
                refseq_id = refseq_id.replace(';','')
                new_isoform.addLocusNames("RefSeq", refseq_id)
                added_gencode_or_refseq = True

            if (added_gencode_or_refseq):
                annotator.setHGNCIfPossible(new_isoform)

            new_isoform.addLocusNames("AceView", gene_id)
            isoform_models[transcript_id] = new_isoform

            if (start_supported == "N"):
                new_isoform.setFivePrimeIncomplete()
            elif (start_supported == "Y"):
                new_isoform.setFivePrimeComplete()
            else:
                assert (start_supported == "U")

            if (end_supported == "N"):
                new_isoform.setThreePrimeIncomplete()
            elif (end_supported == "Y"):
                new_isoform.setThreePrimeComplete()
            else:
                assert (end_supported == "U")

        if (fields[2] == "exon"):
            chromosome = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            strand = fields[6]
            mo = reExonAnnotation.match(fields[-1])
            assert (mo != None)
            gene_id, transcript_id = mo.groups()
                
            try:
                isoform_models[transcript_id].addExon(chromosome, start, stop, strand)
            except AssertionError as ae:
                if (ae.message == "Inconsistent exon strand"):
                    # This happens because of liftover from hg19 to hg38
                    loci_to_remove.add(gene_id)
                elif (ae.message == "Inconsistent chromosome"):
                    pdb.set_trace()
                    print >> sys.stderr, transcript_id, chromosome, start, stop, strand
    ip.close()
    
    transcripts_removed_per_locus = defaultdict(list)
    for transcript_id, isoform in isoform_models.items():
        locus = isoform.getSrcDBLocus()
        if (locus in loci_to_remove):
            transcripts_removed_per_locus[locus].append(transcript_id)
            
    if (len(loci_to_remove) > 0):
        print >> sys.stderr, "WARNING: removing transcripts from the following %d loci due to inconsistent exon strands:" % len(loci_to_remove)
        print >> sys.stderr, "\n".join(map(lambda l: "\t%s\t(%d transcripts)" % (l, len(transcripts_removed_per_locus[l])), loci_to_remove))

    for transcript_id in chain.from_iterable(transcripts_removed_per_locus.values()):
        print >> sys.stderr, "WARNING: removing %s, inconsistent exon strands" % transcript_id
        del isoform_models[transcript_id]

    return isoform_models.values()


def parseGENCODE(db, annotator):
    isoform_models = {}
    ip = openDatabaseFile(db)

    for line in ip:
        fields = line.strip().split("\t")
        #if (fields[0] != "chr8"):
        #    continue

        if (line[0] != '#' and (fields[2] == "transcript" or fields[2] == "exon")):
            chromosome = fields[0]
            strand = fields[6]

            attributes_dict = dict(map(lambda x: tuple(x.split('=')), fields[-1].split(';')))

            #transcript_id = attributes_dict['transcript_id'].split('.')[0]
            
            if (fields[2] == "transcript"):
                transcript_id = attributes_dict['ID'] #.split('.')[0]
                gene_id = attributes_dict['gene_id'].split('.')[0]

                #if (isoform_models.has_key(transcript_id)):
                #    print >> sys.stderr, "WARNING: duplicated transcript ID -> %s" % transcript_id

                try:
                    assert (not isoform_models.has_key(transcript_id))
                except AssertionError as ae:
                    pdb.set_trace()
                    
                new_isoform = RNAIsoform("GENCODE", gene_id, transcript_id, chromosome, strand)
                new_isoform.setFivePrimeComplete()  
                new_isoform.setThreePrimeComplete()

                if (annotator != None):
                    annotator.setHGNCIfPossible(new_isoform)
                if (attributes_dict.has_key("gene_name")):
                    new_isoform.addLocusNames("GENCODE", attributes_dict["gene_name"])
                if (attributes_dict.has_key("tag")):
                    if ("mRNA_start_NF" in attributes_dict["tag"]):
                        new_isoform.setFivePrimeIncomplete()
                    if ("mRNA_end_NF" in attributes_dict["tag"]):
                        new_isoform.setThreePrimeIncomplete()

                isoform_models[transcript_id] = new_isoform
            else:
                transcript_id = attributes_dict['Parent'] #.split('.')[0]
                assert (isoform_models.has_key(transcript_id))
                isoform_models[transcript_id].addExon(chromosome, int(fields[3]), int(fields[4]), strand)

    ip.close()
    return isoform_models.values()


def parseUCSC(db, annotator):
    isoform_models = {}
    ip = openDatabaseFile(db)

    reExonAnnotation = re.compile("gene_id .(\S+).. transcript_id .(\S+)\.\d.. ")
    for line in ip:
        fields = line.strip().split("\t")
        #if (fields[0] != "chr8"):
        #    continue

        if (line[0] != '#' and fields[2] == "exon"):
            chromosome = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            strand = fields[6]
            mo = reExonAnnotation.match(fields[-1])
            assert (mo != None)
            gene_id, transcript_id = mo.groups()
            transcript_id = transcript_id.split('.')[0]
            if (not isoform_models.has_key(transcript_id)):
                assert (not isoform_models.has_key(transcript_id))
                new_isoform = RNAIsoform("UCSC", None, transcript_id, chromosome, strand)
                if (annotator != None):
                    annotator.setHGNCIfPossible(new_isoform)
                # TODO: new_isoform.addLocusNames(database, alt_names_list)
                isoform_models[transcript_id] = new_isoform
            isoform_models[transcript_id].addExon(chromosome, start, stop, strand)

    ip.close()

    return isoform_models.values()


def parseALTSCAN(db, annotator):
    isoform_models = {}
    ip = openDatabaseFile(db)

    for line in ip:
        fields = line.strip().split("\t")

        if (line[0] != '#' and fields[2] == "CDS"):
            chromosome = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            strand = fields[6]
            annot_part = fields[-1].replace(';','')
            annot_part = annot_part.split()
            annot_dict = dict([(annot_part[i],annot_part[i+1]) for i in xrange(0,len(annot_part),2)])
            gene_id = annot_dict["gene_id"]
            transcript_id = annot_dict["transcript_id"]

            if (gene_id.startswith("VHC_")): # Use only the medium confidence transcripts, which include the very high confidence
                continue
            elif (not isoform_models.has_key(transcript_id)):
                new_isoform = RNAIsoform("ALTSCAN", gene_id, transcript_id, chromosome, strand)
                new_isoform.setFivePrimeIncomplete()
                new_isoform.setThreePrimeIncomplete()
                if (annot_dict.has_key("HGNC")):
                    new_isoform.setHGNCID(annot_dict["HGNC"])
                if (annot_dict.has_key("RefSeq")):
                    new_isoform.addLocusNames("RefSeq", annot_dict["RefSeq"])
                isoform_models[transcript_id] = new_isoform
            isoform_models[transcript_id].addExon(chromosome, start, stop, strand)

    ip.close()

    return isoform_models.values()


def parseRefSeq(db, annotator):
    isoform_models = {}
    is_pseudogene = set()
    
    ip = openDatabaseFile(db)

    features_to_parse = set(["primary_transcript", "gene", "transcript", "mRNA", "exon",
                             "miRNA", "lnc_RNA", "snRNA", "snoRNA", "ncRNA", "tRNA", "rRNA", "antisense_RNA",
                             "telomerase_RNA", "vault_RNA", "Y_RNA", "RNase_MRP_RNA", "RNase_P_RNA",
                             "C_gene_segment", "D_gene_segment", "J_gene_segment", "V_gene_segment"])

    mRNAs_to_remove = set()
    
    for line in ip:
        fields = line.strip().split("\t")
        #if (fields[0] != "chr8"):
        #    continue

        if (line[0] != '#' and (fields[2] in features_to_parse)):
            transcript_id, gene_id = None, None
            ID, Parent = None, None  # In the GFF3 attribute field
            chromosome = fields[0]
            source = fields[1]
            strand = fields[6]

            attributes_dict = dict(map(lambda x: tuple(x.split('=')), fields[-1].split(';')))

            if (fields[2] == "exon"):
                Parent = attributes_dict["Parent"]
                try:
                    assert (isoform_models.has_key(Parent))
                except AssertionError as ae:
                    pdb.set_trace()
                    
                isoform_models[Parent].addExon(chromosome, int(fields[3]), int(fields[4]), strand)

                if (source == "Gnomon" and attributes_dict.has_key("exception") and "discrepancy" in attributes_dict["exception"]):
                    mRNAs_to_remove.add(Parent)
            else:
                gene_id = attributes_dict['gene']
                ID = attributes_dict["ID"]
                if (attributes_dict.has_key("pseudo") and attributes_dict["pseudo"] == "true"):
                    transcript_id = gene_id
                elif (not attributes_dict.has_key("transcript_id")):
                    transcript_id = gene_id
                else:
                    transcript_id = attributes_dict["transcript_id"].split('.')[0]
                        
                assert (not isoform_models.has_key(ID))
                new_isoform = RNAIsoform("RefSeq", gene_id, transcript_id, chromosome, strand)
                new_isoform.setFivePrimeComplete()
                new_isoform.setThreePrimeComplete()

                if (attributes_dict.has_key('Dbxref')):
                    db_xref = {}
                    for Db, Db_acc in map(lambda y: tuple(y.split(':', 1)), attributes_dict['Dbxref'].split(',')):
                        if (Db == "HGNC" and not Db_acc.endswith("HGNC")): # Some entries are corrupted with "HGNC:HGNC:HGNC" in the RefSeq file
                            new_isoform.setHGNCID(Db_acc)

                if (annotator != None):
                    annotator.setHGNCIfPossible(new_isoform) # Performs check against data from HGNC, so not superfluous
                if (attributes_dict.has_key("gene")):
                    new_isoform.addLocusNames("RefSeq", attributes_dict["gene"])
                isoform_models[ID] = new_isoform

                # The RefSeq database has some "special" cases that are a problem to handle
                if (source == "Gnomon" and attributes_dict.has_key("exception") and "discrepancy" in attributes_dict["exception"]):
                    mRNAs_to_remove.add(ID)
                if (attributes_dict.has_key("Note") and "frameshift" in attributes_dict["Note"]):
                    mRNAs_to_remove.add(ID)                    
    ip.close()

    # Because some pseudo genes don't consistently obey the gene->mRNA->exon hierarchy (the mRNA is *sometimes* missing),
    # some pseudogene isoforms don't have exons. Remove these cases.
    # Consider using https://pypi.python.org/pypi/gff3.
    for transcript_id, isoform_instance in isoform_models.items():
        if (isoform_instance.getNumberOfExons() == 0):
            mRNAs_to_remove.add(transcript_id)

    for transcript_id in mRNAs_to_remove:
        del isoform_models[transcript_id]

    return isoform_models.values()


def parseSIB(db, annotator):
    assert (annotator != None), "Annotator is null in call to parseSIB()"

    isoform_models = {}
    ip = openDatabaseFile(db)

    reExonAnnotation = re.compile("gene_id .(\S+).. transcript_id .(\S+).. ")
    reGENCODE = re.compile("GENCODE \"(\S+)\";")
    reRefSeq = re.compile("RefSeq \"(\S+)\";")

    for line in ip:
        fields = line.strip().split("\t")

        if (line[0] != '#' and fields[2] == "exon"):
            chromosome = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            strand = fields[6]
            mo = reExonAnnotation.match(fields[-1])
            assert (mo != None)
            gene_id, transcript_id = mo.groups()
            assert (gene_id == transcript_id)
            gene_id = gene_id.rsplit('.',1)[0]
            if (not isoform_models.has_key(transcript_id)):
                new_isoform = RNAIsoform("SIB", gene_id, transcript_id, chromosome, strand)

                added_gencode_or_refseq = False
                mo = reGENCODE.search(fields[-1])
                if (mo != None):
                    gencode_id = mo.group(1)
                    new_isoform.addLocusNames("GENCODE", gencode_id)
                    added_gencode_or_refseq = True

                mo = reRefSeq.search(fields[-1])
                if (mo != None):
                    refseq_id = mo.group(1)
                    new_isoform.addLocusNames("RefSeq", refseq_id)
                    added_gencode_or_refseq = True

                if (added_gencode_or_refseq):
                    annotator.setHGNCIfPossible(new_isoform)

                isoform_models[transcript_id] = new_isoform
            
            isoform_models[transcript_id].addExon(chromosome, start, stop, strand)

    ip.close()

    return isoform_models.values()


def parseHInv(db, annotator):
    isoform_models = {}
    ip = openDatabaseFile(db)

    #remRNAAnnotation = re.compile("ID=(\S+).Name.*Note=(\S*)\/.+Parent=(HIX\d+)")
    #reExonAnnotation = re.compile("Parent=(\S+)\s*$")

    for line in ip:
        fields = line.strip().split("\t")

        if (line[0] != '#' and fields[2] == "mRNA"):
            try:
                attributes_dict = dict(map(lambda x: tuple(x.split('=')), filter(lambda y:'=' in y, fields[-1].split(';'))))
            except ValueError as ve:
                print ve.message
                pdb.set_trace()

            transcript_id = attributes_dict["ID"]
            gene_name = attributes_dict["Note"].split('/')[0]
            parent_id = attributes_dict["Parent"]

            assert (not isoform_models.has_key(transcript_id))

            new_isoform = RNAIsoform("HInv", parent_id, transcript_id, fields[0], fields[6])

            added_gencode_or_refseq = False
            if (attributes_dict.has_key("GENCODE")):
                gencode_id = attributes_dict["GENCODE"]
                new_isoform.addLocusNames("GENCODE", gencode_id)
                added_gencode_or_refseq = True

            if (attributes_dict.has_key("RefSeq")):
                refseq_id = attributes_dict["RefSeq"]
                new_isoform.addLocusNames("RefSeq", refseq_id)
                added_gencode_or_refseq = True

            if (added_gencode_or_refseq):
                annotator.setHGNCIfPossible(new_isoform)

            if (gene_name != ''):
                new_isoform.addLocusNames("HInv", gene_name)
            isoform_models[transcript_id] = new_isoform

        elif (line[0] != '#' and fields[2] == "exon"):
            chromosome = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            strand = fields[6]
            attributes_dict = dict(map(lambda x: tuple(x.split('=')), filter(lambda y:'=' in y, fields[-1].split(';'))))
            transcript_id = attributes_dict["Parent"]
            isoform_models[transcript_id].addExon(chromosome, start, stop, strand)
                
    ip.close()

    return isoform_models.values()


def parseRNAcentral(db): #, annotator
    num_other_db_skipped = 0
    isoform_models = {}
    ip = openDatabaseFile(db)

    #reExonAnnotation = re.compile("ID\s+.(\S+)..Name\s+.(\S+).$")
    for line in ip:
        fields = line.strip().split("\t")
        #if (fields[0] != "chr8"):
        #    continue

        if (line[0] != '#'):
            chromosome = fields[0]
            strand = fields[6]
            attributes_dict = dict(map(lambda x: tuple(x.split()), fields[-1].split(';')))

            if (fields[2] == "transcript"):
                #mo = reExonAnnotation.match(fields[-1])
                #assert (mo != None)
                #transcript_id, gene_id = mo.groups()
                transcript_id = attributes_dict["ID"].strip("\"")
                if (transcript_id.startswith("NR_") or transcript_id.startswith("ENST0")):
                    num_other_db_skipped += 1
                    continue # This transcript will be added when reading RefSeq and/or GENCODE databases
                gene_id = attributes_dict["Name"].strip("\"")
                assert (not isoform_models.has_key(transcript_id))
                new_isoform = RNAIsoform("RNAcentral", gene_id, transcript_id, chromosome, strand)
                #annotator.setHGNCIfPossible(new_isoform)
                # TODO: new_isoform.addLocusNames(database, alt_names_list)
                isoform_models[transcript_id] = new_isoform
            else:
                assert(fields[2] == "noncoding_exon")
                transcript_id = attributes_dict["Parent"].strip("\"")
                if (isoform_models.has_key(transcript_id)): # It wouldn't if RNA is from RefSeq or GENCODE/ENSEMBL
                    try:
                        isoform_models[transcript_id].addExon(chromosome, int(fields[3]), int(fields[4]), strand)
                    except KeyError as ke:
                        print >> sys.stderr, "WARNING: no isoform set for %s. Repairing." % transcript_id
                        gene_id = attributes_dict["Name"].strip("\"")
                        isoform_models[transcript_id] = RNAIsoform("RNAcentral", gene_id, transcript_id, chromosome, strand)
                        isoform_models[transcript_id].addExon(chromosome, int(fields[3]), int(fields[4]), strand)
    ip.close()

    print >> sys.stderr, "INFO: skipped %d isoforms that derive from RefSeq or GENCODE" % num_other_db_skipped
    return isoform_models.values()


def parseCGDB(db):
    isoform_models = {}
    ip = openDatabaseFile(db)

    reExonAnnotation = re.compile("gene_id .(\S+).. transcript_id .(.+)\";$")
    for line in ip:
        fields = line.strip().split("\t")
        #if (fields[0] != "chr8"):
        #    continue

        if (line[0] != '#' and fields[2] == "exon"):
            chromosome = fields[0]
            start = int(fields[3])
            stop = int(fields[4])
            strand = fields[6]
            mo = reExonAnnotation.match(fields[-1])
            assert (mo != None)
            gene_id, transcript_id = mo.groups()
            if (not isoform_models.has_key(transcript_id)):
                new_isoform = RNAIsoform("hESC", gene_id, transcript_id, chromosome, strand)
                # TODO: new_isoform.addLocusNames(database, alt_names_list)
                isoform_models[transcript_id] = new_isoform

            isoform_models[transcript_id].addExon(chromosome, start, stop, strand)

    ip.close()

    return isoform_models.values()


def parseFANTOM(db):
    TSSs = []
    ip = openDatabaseFile(db)
    
    if ("permissive" in db):
        source_label = "permissive"
    elif ("robust" in db):
        source_label = "robust"
    elif ("phase1and2combined" in db):
        source_label = "phase1and2combined"
    else:
        print >> sys.stderr, "ERROR: unrecognized FANTOM5 source file -> %s. Exiting." % db
        sys.exit(1)
        
    # BED format:
    #chr8    629191  629220  chr8:564571..564600,+   2398    +
    for line in ip:
        fields = line.strip().split("\t")
        #if (fields[0] != "chr8"):
        #    continue
        chromosome = fields[0]
        start = int(fields[1]) + 1
        stop = int(fields[2])
        tss_label = fields[3]   # FANTOM5 data had to be lifted over to hg38, so the label is based on hg19 coordinates. And, is possibly concat'ed ids.
        strand = fields[5]
        TSSs.append( TSS("FANTOM", source_label, tss_label, chromosome, start, stop, strand) )

    ip.close()
    return TSSs


def parseAPADB(db):
    polyA_sites = []
    ip = openDatabaseFile(db)
    
    # BED format:
    #chr8   198356  198358  RPH3AL.2        8       -       Intron  198358.0        198358  miR-149,miR-124/124ab/506,miR-141/200a
    for line in ip:
        fields = line.split("\t")
        #if (fields[0] != "chr8"):
        #    continue
        chromosome = fields[0]
        start = int(fields[1]) + 1
        stop = int(fields[2])
        pAs_label = fields[3]   # FANTOM5 data had to be lifted over to hg38, so the label is based on hg19 coordinates
        strand = fields[5]
        cluster_median = int(round(float(fields[7])))
        cluster_mode = int(fields[8])
        polyA_sites.append( PolyA("APADB", pAs_label, chromosome, start, stop, strand, cluster_median, cluster_mode) )

    ip.close()
    return polyA_sites
