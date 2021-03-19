#!/usr/bin/env python

import itertools
import sys

if len(sys.argv) < 5:
	print "usage: %s <dosefile> <idfile> <pedfile> <mapfile>" % sys.argv[0]
	sys.exit(2)

sys.stderr.write("dose file: %s\n" % sys.argv[1])
sys.stderr.write("id file: %s\n" % sys.argv[2])
with open(sys.argv[1],'rU') as doseFile:
	with open(sys.argv[2],'rU') as idFile:
		sys.stderr.write("writing map file: %s\n" % sys.argv[4])
		n = 0
		with open(sys.argv[4],'wb') as mapFile:
			for snp in doseFile.next().strip().split():
				n += 1
				if n > 1: # skip PHENOTYPE column
					mapFile.write("0 %s 0\n" % snp)
		#with mapFile
		sys.stderr.write("... OK: %d SNPs\n" % n)
		
		sys.stderr.write("writing ped file: %s\n" % sys.argv[3])
		n = 1 # the header line we got from .next() above
		idHeader = idFile.next().strip()
		if idHeader != "ID":
			sys.exit("ERROR: unknown id file header; expected 'ID', found: %s\n" % idHeader)
		with open(sys.argv[3],'wb') as pedFile:
			for line in doseFile:
				n += 1
				inCols = line.strip().split()
				outCols = (('0' if f < 0.5 else ('2' if f > 1.5 else '1')) for f in (float(c) for c in inCols[1:]))
				try:
					iid = idFile.next().strip()
					prefix = [iid,iid,'0','0','0',inCols[0]]
				except StopIteration:
					sys.exit("ERROR: id file ends at line #%d, dose file continues\n" % n)
				pedFile.write(" ".join(itertools.chain(prefix, outCols))+"\n")
		sys.stderr.write("... OK: %d lines\n" % (n,))
		
		try:
			idFile.next()
			sys.stderr.write("WARNING: dose file ends at line #%d, id file continues\n" % n)
		except StopIteration:
			pass
	#with idFile
#with doseFile
