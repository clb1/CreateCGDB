#!/usr/bin/env python

from collections import defaultdict
import re
import gzip
import sys


# mRNA  AceView Gene    Entrez Gene ID (if any)         included RefSeq (if any) 
def parseAliases(aliases_tsv):
    transcript_aliases = defaultdict(set)
    ip = gzip.open(aliases_tsv, 'rb')
    for line in ip:
        line = line.strip()
        entrezID, refseqID = "",""
        if (line != "" and line[0] != '#' and line[0] != '/' ):
            transcript_id, gene_id, entrez, refseq = line.replace("\"","").split("\t")
            if (entrez != "NULL"):
                if (";" in entrez):
                    mult_entrez = entrez.replace("\\;","").split()
                    transcript_aliases[transcript_id].add(",".join(map(lambda x: "Entrez:%s" % x, mult_entrez)))
                else:
                    transcript_aliases[transcript_id].add("Entrez:%s" % entrez)
            if (refseq != "NULL"):
                transcript_aliases[transcript_id].add("RefSeq:%s" % refseq)
    ip.close()

    return transcript_aliases

    
def modifyAndWrite(transcript_aliases, output_gff3):
    if (output_gff3.endswith(".gz")):
        op = gzip.open(output_gff3, 'wb')
    else:
        op = open(output_gff3, 'w')
    
    reGeneAttributes = re.compile("ID=(\S+);Name")
    reTranscriptAttributes = re.compile("ID=(\S+);Parent")
    
    for line in sys.stdin:
        if (line[0] == '#'):
            op.write(line)
        else:
            fields = line.strip().split("\t")
            fields[-1] = fields[-1].replace(";Gene_type=cDNA_supported","")

            #if (fields[2] == "gene"):
            #    mo = reGeneAttributes.match(fields[-1])
            #    assert (mo != None)
            #    gene_name = mo.group(1)
            #    if (gene_aliases.has_key(gene_name)):
            #        fields[-1] = "%s;Alias=%s" % (fields[-1], gene_aliases[gene_name])
            #    fields[0] = "###\n%s" % fields[0]
            if (fields[2] == "transcript"):
                mo = reTranscriptAttributes.match(fields[-1])
                assert (mo != None)
                transcript_name = mo.group(1)
                if (transcript_aliases.has_key(transcript_name)):
                    fields[-1] = "%s;Alias=%s" % (fields[-1], ",".join(transcript_aliases[transcript_name]))

            op.write("%s\n" % "\t".join(fields))
    op.close()


if (__name__ == "__main__"):
    aliases_tsv, output_gff3 = sys.argv[1:]

    transcript_aliases = parseAliases(aliases_tsv)
    modifyAndWrite(transcript_aliases, output_gff3)
    
    sys.exit(0)
    
