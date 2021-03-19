#!/usr/bin/env python

import sys
import collections


if len(sys.argv) < 4:
	print "usage: %s <IBC-file> <models-file> <output-file>" % sys.argv[0]
	sys.exit(2)

print "loading IBC data from '%s' ..." % sys.argv[1]
locusSNP = collections.defaultdict(lambda: return "?")
with open(sys.argv[1],'rU') as ibcFile:
	for line in ibcFile:
		words = line.strip().lower().split()
		locus = (words[1],words[2])
		if locus in locusSNP:
			if locusSNP[locus] != words[0]:
				print "WARNING: ambiguous SNP/position mapping, only the first will be used:"
				print "%s\t%s\t%s" % (locusSNP[locus],words[1],words[2])
				print "%s\t%s\t%s" % (words[0],words[1],words[2])
		else:
			locusSNP[locus] = words[0]
	#foreach line
#with ibcFile
print "... OK: %d SNP/position mappings" % len(locusSNP)

print "converting models file '%s' ..." % sys.argv[2]
with open(sys.argv[2],'rU') as modelFile:
	line = modelFile.next().strip()
	if line != "chr1	pos1	chr2	pos2	score":
		print "ERROR: unexpected model file header: %s" % line
		sys.exit(1)
	with open(sys.argv[3],'w') as outputFile:
		outputFile.write("snp1	snp2	score\n")
		n = 0
		for line in modelFile:
			n += 1
			words = line.strip().lower().split()
			locus1 = (words[0],words[1])
			locus2 = (words[2],words[3])
			outputFile.write("%s\t%s\t%s\n" % (locusSNP[locus1],locusSNP[locus2],words[4]))
		#foreach line
	#with outputFile
#with modelFile
print "... OK: %d models converted and written to '%s'" % (n,sys.argv[3])

