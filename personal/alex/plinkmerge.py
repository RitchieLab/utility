#!/usr/bin/env python

import sys
import os.path

if len(sys.argv) < 3:
	sys.stderr.write(
"""Usage: %s <file1> <file2>

genotype and phenotype files may be specified in either order
""" % os.path.basename(sys.argv[0]))
	sys.exit(1)

with open(sys.argv[1], 'r') as file1:
	with open(sys.argv[2], 'r') as file2:
		gheader = "FID1\tIID1\tFID2\tIID2\tRT\tEZ\tZ0\tZ1\tZ2\tPI_HAT\tPHE\tDST\tPPC\tRATIO"
		pheader = "FID\tIID\tCENTER\tSTUDY\tRH_CASE\tDECADE_BIRTH"
		
		header1 = "\t".join(file1.next().split())
		header2 = "\t".join(file2.next().split())
		
		if header1.startswith(gheader) and header2.startswith(pheader):
			sys.stderr.write("identified genotype and phenotype files\n")
			gfile = file1
			pfile = file2
		elif header1.startswith(pheader) and header2.startswith(gheader):
			sys.stderr.write("identified phenotype and genotype files\n")
			gfile = file2
			pfile = file1
		else:
			sys.stderr.write("unrecognized file header(s)")
			sys.exit(1)
		
		sys.stderr.write("processing phenotype file ...\n")
		center = {}
		decade_birth = {}
		study = {}
		n = 0
		for line in pfile:
			tok = line.split()
			if len(tok) >= 6:
				n += 1
				iid = tok[1]
				center[iid] = tok[2]
				decade_birth[iid] = tok[5]
				study[iid] = tok[3]
		sys.stderr.write("... complete: %d lines\n" % n)
		
		sys.stderr.write("processing genotype file ...\n")
		sys.stdout.write("FID1\tIID1\tCENTER1\tDECADE_BIRTH1\tSTUDY1\tFID2\tIID2\tCENTER2\tDECADE_BIRTH2\tSTUDY2\tRT\tEZ\tZ0\tZ1\tZ2\tPI_HAT\tPHE\tDST\tPPC\tRATIO\n")
		n = 0
		for line in gfile:
			tok = line.split()
			if len(tok) >= 5:
				n += 1
				iid1 = tok[1]
				iid2 = tok[3]
				out = tok[0:2]
				out.extend([center[iid1], decade_birth[iid1], study[iid1]])
				out.extend(tok[2:4])
				out.extend([center[iid2], decade_birth[iid2], study[iid2]])
				out.extend(tok[4:])
				sys.stdout.write("\t".join(out))
				sys.stdout.write("\n")
		sys.stderr.write("...  complete: %d lines\n" % n)
	#with file2
#with file1
