#!/usr/bin/env python

from collections import defaultdict
import re
import sys

if (__name__ == "__main__"):
	failed_elements = defaultdict(set)
	
	with open(sys.argv[1],'r') as ip:
		for line in ip:
			element, ID = line.strip().split("\t")
			failed_elements[element].add(ID)
			failed_elements["ALL"].add(ID)
			
	with open(sys.argv[2], 'r') as ip:
		referenced_genes = set(ip.read().split('\n'))
		
	reParent = re.compile("Parent=(\S+)$")
	num_elements_filtered = 0
	for line in sys.stdin:
		if (line[0] == '#'):
			continue
		line = line.strip()
		line_elems = line.split("\t")
		element_type = line_elems[2]
		mo = reParent.match(line_elems[-1])
		if (mo != None):
			parentID = mo.group(1)
		else:
			parentID = None
		ID = line_elems[8].split(';')[0][3:]
		chrom_ID = "%s\t%s" % (line_elems[0], ID)
		if (not ID in failed_elements[element_type] and not parentID in failed_elements["ALL"]):
			if (element_type=="gene" and not chrom_ID in referenced_genes):
				num_elements_filtered += 1	
			else:
				print >> sys.stdout, line	
		else:
			num_elements_filtered += 1

	print >> sys.stderr, "INFO: filtered %d elements" % num_elements_filtered
	sys.exit(0)
	
