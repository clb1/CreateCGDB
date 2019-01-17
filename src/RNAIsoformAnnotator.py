import gzip
import sys
from RNAIsoform_python2 import RNAIsoform


class RNAIsoformAnnotator(object):
    def __init__(self, hgnc_data_tsv = None):
        self.hgnc_id_to_symbol = {}
        self.to_hgnc_id = {}
        self.all_hgnc_data = {}
        #self.ncbi_NM_and_symbol_to_hgnc = None
        self.ucsc_to_hgnc_ensembl_entrez = None
        self.rnacentral_to_hgnc_ensembl_entrez = None
        #self.aceview_info = {}

        if (hgnc_data_tsv != None):
            self.readHGNCData(hgnc_data_tsv)

        
    # Format: HGNC ID, Approved Symbol,    Approved Name, Status, Locus Type, Locus Group,    Previous Symbols,
    #         Previous Name, Synonyms, Name Synonyms, Chromosome, Accession Numbers, Entrez Gene ID,
    #         Ensembl Gene ID, RefSeq IDs, Entrez Gene ID(supplied by NCBI), RefSeq(supplied by NCBI),
    #         Ensembl ID(supplied by Ensembl), UCSC ID(supplied by UCSC)
    def readHGNCData(self, hgnc_data_tsv):
        hgnc_id_to_symbol = {}
        to_hgnc_id = {}
        with gzip.open(hgnc_data_tsv, 'r') as ip:
            header = ip.readline()
            for line in ip:
                if (not "withdrawn" in line):
                    fields = line.strip("\n").split("\t")
                    standardized_hgnc_id = fields[0].strip()
                    self.all_hgnc_data[standardized_hgnc_id] = fields[1:]

                    hgnc_symbol = fields[1].strip()
                    hgnc_id_to_symbol[standardized_hgnc_id] = hgnc_symbol
                
                    # Map Entrez GeneID to HGNC ID if possible
                    for entrez in [fields[12].strip(), fields[15].strip()]:
                        if (entrez != ''):
                            standardized_entrez_id = "Entrez:%s" % entrez
                            to_hgnc_id[standardized_entrez_id] = standardized_hgnc_id

                    # Map Ensembl Gene to HGNC ID if possible
                    for ensembl in [fields[13].strip(), fields[17].strip()]:
                        if (ensembl != ''):
                            to_hgnc_id[ensembl] = standardized_hgnc_id

                    # Map RefSeq to HGNC ID if possible
                    for refseq in [fields[14].strip().split('.')[0], fields[16].strip().split('.')[0]]:
                        if (refseq != ''):
                            to_hgnc_id[refseq] = standardized_hgnc_id

                    hgnc_ucsc = fields[18].strip()
                    if (hgnc_ucsc != ""):
                        hgnc_ucsc = hgnc_ucsc.split('.')[0]
                        to_hgnc_id[hgnc_ucsc] = standardized_hgnc_id

        self.hgnc_id_to_symbol = hgnc_id_to_symbol
        self.to_hgnc_id = to_hgnc_id


    def readBiomartUCSC(self, biomart_ucsc_tsv):
        biomart_ucsc = defaultdict(set)
        with gzip.open(biomart_ucsc_tsv, 'rb') as ip:
            header = ip.readline()
            for line in ip:
                ensembl_gene, ensembl_transcript, hgnc_id, entrez_id, ucsc_id = line.strip("\n").split("\t")
                if (ucsc_id != ""):
                    ucsc_id = ucsc_id.split('.')[0]
                    biomart_ucsc[ucsc_id].add( (ensembl_gene, ensembl_transcript, hgnc_id, entrez_id) )
        self.ucsc_to_hgnc_ensembl_entrez = biomart_ucsc


    def readBiomartRNACentral(self, biomart_RNACentral_tsv):
        biomart_rnacentral = defaultdict(set)
        with gzip.open(biomart_RNACentral_tsv, 'rb') as ip:
            header = ip.readline()
            for line in ip:
                ensembl_gene, ensembl_transcript, hgnc_id, rnacentral_id, entrez_id = line.strip("\n").split("\t")
                if (rnacentral_id != ""):
                    biomart_rnacentral[rnacentral_id].add( (ensembl_gene, ensembl_transcript, hgnc_id, entrez_id) )
        self.rnacentral_to_hgnc_ensembl_entrez = biomart_rnacentral


    # Format:
    # mRNA  AceView Gene    Entrez Gene ID (if any)         included RefSeq (if any) 
    def readAceViewAliases(self, aliases_tsv):
        ip = gzip.open(aliases_tsv, 'rb')
        for line in ip:
            if (not line[0] in "\n#/"):
                line = line.replace("\"","")
                transcript_id, gene_id, entrez, refseq = line.strip("\n").split("\t")
                if (";" in entrez): # Transcript spans multiple Entrez genes
                    pass
                    #standardized_entrez_id = set(map(lambda x: "Entrez:%s" % x, entrez.split("\; ")))
                    #conflicting_hgnc_ids = set()
                    #for entrez in standardized_entrez_id:
                    #    if (to_hgnc_id.has_key(entrez)):
                    #        hgnc_id = to_hgnc_id[entrez]
                    #        conflicting_hgnc_ids.add(hgnc_id)
                    #if (len(conflicting_hgnc_ids) > 1):
                    #    hgnc_id = None
                    #elif (len(conflicting_hgnc_ids) == 1):
                    #    hgnc_id = list(conflicting_hgnc_ids)[0]
                    #    conflicting_hgnc_ids = None
                    #else:
                    #    hgnc_id = None
                    #    conflicting_hgnc_ids = None
                    #refseq = "" if (refseq != "NULL") else refseq
                else:
                    hgnc_ids = set()
                    if (refseq != "NULL"):
                        if (self.to_hgnc_id.has_key(refseq)):
                            hgnc_ids.add( self.to_hgnc_id[refseq] )
                            
                    elif (entrez != "NULL"):
                        standardized_entrez_id = "Entrez:%s" % entrez
                        if (self.to_hgnc_id.has_key(standardized_entrez_id)):
                            hgnc_ids.add( self.to_hgnc_id[standardized_entrez_id] )

                    if (len(hgnc_ids) == 1):
                        self.to_hgnc_id[transcript_id] = list(hgnc_ids)[0]

        ip.close()


    def setHGNCIfPossible(self, isoform):
        src_db = isoform.getSrcDB()

        if (src_db == "GENCODE"):
            src_db_locus = isoform.getSrcDBLocus()
            if (self.to_hgnc_id.has_key(src_db_locus)):
                hgnc_id = self.to_hgnc_id[src_db_locus]
                hgnc_symbol = self.hgnc_id_to_symbol[hgnc_id]
                isoform.setHGNC(hgnc_id, hgnc_symbol)

        elif (src_db in ["AceView", "UCSC"]):
            transcript_id = isoform.getTranscriptIDInSrcDB()
            if (self.to_hgnc_id.has_key(transcript_id)):
                hgnc_id = self.to_hgnc_id[transcript_id]
                hgnc_symbol = self.hgnc_id_to_symbol[hgnc_id]
                isoform.setHGNC(hgnc_id, hgnc_symbol)

        elif (src_db == "RefSeq"):
            updated = False
            hgnc_id = isoform.getHGNCID()
            hgnc_symbol = isoform.getHGNCSymbol()

            if (hgnc_id == None):
                transcript_id = isoform.getTranscriptIDInSrcDB()
                if (self.to_hgnc_id.has_key(transcript_id)):
                    hgnc_id = self.to_hgnc_id[transcript_id]
                    updated = True

            if (hgnc_id != None and hgnc_symbol == None and self.hgnc_id_to_symbol.has_key(hgnc_id)):
                hgnc_symbol = self.hgnc_id_to_symbol[hgnc_id]
                updated = True

            if (updated):
                isoform.updateHGNC(hgnc_id, hgnc_symbol)


        # Some H-Inv,SIB,Aceview database isoforms have custom GENCODE/RefSeq locus annotations, which can be used to get HGNC
        if (src_db == "HInv" or src_db == "SIB" or src_db == "AceView"): 
            gencode_locus = isoform.getLocusSymbolFromDatabase("GENCODE")
            refseq_locus = isoform.getLocusSymbolFromDatabase("RefSeq")
            if (gencode_locus != None):
                gencode_locus = gencode_locus.split('.')[0]
                if (self.to_hgnc_id.has_key(gencode_locus)):
                    hgnc_id = self.to_hgnc_id[gencode_locus]
                    hgnc_symbol = self.hgnc_id_to_symbol[hgnc_id]
                    isoform.setHGNC(hgnc_id, hgnc_symbol)
            elif (refseq_locus != None and self.to_hgnc_id.has_key(refseq_locus)):
                hgnc_id = self.to_hgnc_id[gencode_locus]
                hgnc_symbol = self.hgnc_id_to_symbol[hgnc_id]
                isoform.setHGNC(hgnc_id, hgnc_symbol)



    def resolveIntersectResults(self, isoform, isoform_to_group_map, intersect_results):
        if (len(intersect_results) == 1 and intersect_results[0][2] == None): # Isoform overlaps no other
            assert (not isoform_to_group_map.has_key(isoform1))
            ig = IsoformGroup()
            ig.addIsoform(isoform1)
            isoform_to_group_map[isoform1] = ig
        else:
            pass
    
