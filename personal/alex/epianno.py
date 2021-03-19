#!/usr/bin/env python

import apsw
import os
import sys


if len(sys.argv) <= 2:
	print "usage: %s <dbfile> <mapfile>" % (sys.argv[0],)
	sys.exit(2)

dbfilename = sys.argv[1]
db = apsw.Connection(dbfilename)
cursor = db.cursor()
mapfilename = sys.argv[2]
#threshold = float(sys.argv[3])


def getChromData(chr):
	raise Exception("unused")
	chunkData = list()
	runningState = 0
	sql = "SELECT *,(chunk-1)*200,chunk*200-1,1 FROM epigenome WHERE chr=? ORDER BY chunk"
	for row in cursor.execute(sql, (chr,)):
		if row[2+row[2]] >= threshold:
			if runningState == row[2]:
				for c in xrange(3,28):
					chunkData[-1][c] += row[c]
				for c in xrange(28,53):
					chunkData[-1][c] += row[c]
				chunkData[-1][-3] = min(chunkData[-1][-3], row[-3])
				chunkData[-1][-2] = max(chunkData[-1][-2], row[-2])
				chunkData[-1][-1] += 1
				chunkData.append(chunkData[-1])
			else:
				runningState = row[2]
				chunkData.append(list(row))
		else:
			runningState = 0
			chunkData.append(list(row))
	#for row
	last = None
	for r,row in enumerate(chunkData):
		if last != row[1]:
			last = row[1]
			if row[-1] > 1:
				for c in xrange(3,28):
					chunkData[r][c] /= float(row[-1])
				for c in xrange(28,53):
					chunkData[r][c] /= float(row[-1])
				#for prob col
			#if merged
		#skip repeats
	#for row
	return chunkData
#getChromData()


with open(mapfilename,'rU') as mapfile:
	states = range(1,16)
	stateTally = list(0 for s in states)
	sql = "SELECT %s FROM epigenome WHERE chr=? AND chunk=(? / 200)+1;" % (",".join(("s%s" % s) for s in states),)
	
	header = next(mapfile).split()
	header.extend(["START","STOP"])
	header.extend(list(("S%d" % s) for s in states))
	print "\t".join(header)
	
	for line in mapfile:
		line = line.split()
		chm = line[1]
		pos = long(line[2])
		for s in states:
			stateTally[s-1] = 0
		for row in cursor.execute(sql, (chm,pos)):
			stateTally[max(states, key=lambda s: row[s-1])-1] += 1
		print "%s\t%d\t%d\t%s" % (
			"\t".join(line),
			int(pos/200),
			int(pos/200)+199,
			"\t".join(stateTally),
		)
#with mapfile
