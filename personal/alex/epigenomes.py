#!/usr/bin/env python

import apsw
import os
import sys
import zlib


class zopen(object):
	
	def __init__(self, fileName, splitChar="\n", chunkSize=16*1024):
		try:
			self._filePtr = open(fileName,'rb')
		except:
			self._filePtr = None
			raise
		self._splitChar = splitChar
		self._chunkSize = chunkSize
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		self._text = ""
		self._lines = list()
	#__init__()
	
	
	def __del__(self):
		if self._filePtr:
			self._filePtr.close()
			self._filePtr = None
	#__del__()
	
	
	def __enter__(self):
		return self
	#__enter__()
	
	
	def __exit__(self, excType, excVal, excTrace):
		pass
	#__exit__()
	
	
	def __iter__(self):
		return self
	#__iter__()
	
	
	def __next__(self):
		# if lines are still cached from the last read, pop one
		if len(self._lines) > 0:
			return self._lines.pop()
		# if there's data left in the source file, read and decompress another chunk
		if self._dc:
			data = self._dc.unused_data
			if data:
				self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
			elif not self._filePtr:
				raise Exception("cannot read a closed file")
			else:
				data = self._filePtr.read(self._chunkSize)
			if data:
				self._text += self._dc.decompress(data)
				data = None
			else:
				self._text += self._dc.flush()
				self._dc = None
		# if there's no text left, we're done
		if not self._text:
			raise StopIteration
		# split the text into lines
		self._lines = self._text.split(self._splitChar)
		self._text = ""
		# if there's more than one line, store the last to combine with the next chunk
		# (but if there's only one line, and more to read, then keep reading until we get a linebreak)
		if len(self._lines) > 1:
			self._text = self._lines.pop()
		elif self._dc:
			self._text = self._lines.pop()
			self._chunkSize *= 2
			return self.__next__()
		# reverse the remaining lines into a stack and pop one to return
		self._lines.reverse()
		return self._lines.pop()
	#__next__()
	
	
	def next(self):
		return self.__next__()
	#next()
	
	
	def seek(self, offset, whence = 0):
		if not self._filePtr:
			raise Exception("cannot seek a closed file")
		if offset != 0:
			raise Exception("zfile.seek() does not support offsets != 0")
		self._filePtr.seek(0, whence)
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		self._text = ""
		self._lines = list()
	#seek()
	
	
	def close(self):
		if self._filePtr:
			self._filePtr.close()
			self._filePtr = None
	#close()
	
#zopen


if len(sys.argv) < 3:
	print "usage: %s <dbfile> <datapath> [numstates=15]" % (sys.argv[0],)
	sys.exit(2)

dbfilename = sys.argv[1]
pathname = sys.argv[2]
numstates = int(sys.argv[3]) if (len(sys.argv) > 3) else 15
print "opening database '%s' ..." % (dbfilename,)
db = apsw.Connection(dbfilename)
cursor = db.cursor()
cursor.execute("PRAGMA journal_mode = OFF")
cursor.execute("PRAGMA synchronous = OFF")
cursor.execute("""
CREATE TABLE IF NOT EXISTS setting (
	setting VARCHAR(32) PRIMARY KEY NOT NULL,
	value VARCHAR(128) NOT NULL
)
""")
cursor.execute("""
INSERT OR IGNORE INTO setting (setting,value) VALUES
	('numstates',%d)
""" % (numstates,))
cursor.execute("""
CREATE TABLE IF NOT EXISTS tissue (
	tissue_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
	tissue VARCHAR(32) UNIQUE NOT NULL,
	description VARCHAR(128) NOT NULL
)
""")
cursor.execute("""
INSERT OR IGNORE INTO tissue (tissue_id,tissue,description) VALUES
	(1,'lungcarcinoma','A549 EtOH 0.02% Lung Carcinoma'),
	(2,'adrenal','Adrenal'),
	(3,'blood','Blood'),
	(4,'bone','Bone'),
	(5,'brain','Brain'),
	(6,'breast','Breast'),
	(7,'cervix','Cervix'),
	(8,'leukemia','Dnd41 TCell Leukemia'),
	(9,'esc','ESC/ESC-Derived'),
	(10,'fatskin','Fat/Skin'),
	(11,'gi','GI'),
	(12,'gismooth','GI Smooth Muscle'),
	(13,'heart','Heart'),
	(14,'hepcarcinoma','HepG2 Hepatocellular Carcinoma'),
	(15,'imr90cell','IMR90 Cell Line'),
	(16,'ipsc','IPSC'),
	(17,'kidney','Kidney'),
	(18,'liver','Liver'),
	(19,'lung','Lung'),
	(20,'muscle','Muscle'),
	(21,'ovary','Ovary'),
	(22,'pancreas','Pancreas'),
	(23,'placenta','Placenta'),
	(24,'spleen','Spleen'),
	(25,'stromalstem','Stromal Connective Stem cells'),
	(26,'thymus','Thymus'),
	(27,'vascular','Vascular')
""")
cursor.execute("""
CREATE TABLE IF NOT EXISTS epigenome (
	chr TINYINT NOT NULL,
	chunk BIGINT NOT NULL,
	tissue_id INTEGER NOT NULL,
	%s
	PRIMARY KEY (chr,chunk,tissue_id)
)
""" % ("\n\t".join(("s%d FLOAT UNSIGNED NOT NULL," % (s+1,)) for s in xrange(numstates)),))

tissueEIDs = {
	1: [114],
	2: [80],
	3: [62,34,45,33,44,43,39,41,42,40,37,48,38,47,124,116,123,29,31,35,51,50,36,32,46,30],
	4: [129],
	5: [71,74,68,69,72,67,73,70,82,81,125,54,53],
	6: [119,27,28],
	7: [117],
	8: [115],
	9: [2,8,1,15,14,16,3,24,7,9,10,13,12,11,4,5,6],
	10: [63,25,23,126,127,55,56,59,61,57,58],
	11: [106,75,77,79,85,84,109,101,102,92,110,94],
	12: [76,78,103,111],
	13: [83,104,95,105],
	14: [118],
	15: [17],
	16: [20,19,18,21,22],
	17: [86],
	18: [66],
	19: [128,88,96],
	20: [120,121,100,108,107,89,52,90],
	21: [97],
	22: [87,98],
	23: [99,91],
	24: [113],
	25: [26,49],
	26: [112,93],
	27: [122,65],
}

print "path = %s" % (pathname,)
print "numstates = %d" % (numstates,)

for chm in xrange(1,24):
	print "processing chromosome %s ..." % (('X' if (chm==23) else str(chm)),)
	cursor.execute("BEGIN TRANSACTION")
	
	# open all data files
	epiFile = dict()
	for t in tissueEIDs:
		for e in tissueEIDs[t]:
			filename = os.path.join(pathname, "E%03d_%d_coreMarks_chr%s_posterior.txt.gz" % (e, numstates, (('X' if (chm==23) else str(chm)))))
			fileobj = zopen(filename)
			epiFile[e] = fileobj
		#for e
	#for t
	print "  found %d epigenome files" % (len(epiFile),)
	
	# verify file headers
	header = ("E%03d\tchr%s" % (e, (('X' if (chm==23) else str(chm)))))
	subheader = "\t".join(("E%d" % (s+1,)) for s in xrange(numstates))
	for e,fileobj in epiFile.iteritems():
		line = next(fileobj).strip()
		if line != header:
			exit("  ERROR: file #%03d header mismatch" % (epi,))
		line = next(fileobj).strip()
		if line != subheader:
			exit("  ERROR: file #%03d subheader mismatch" % (e,))
	#for e,fileobj
	
	# scan all data files
	chunk = 0
	stateCol = range(3+0,3+numstates)
	stateProbs = [list() for s in xrange(numstates)]
	tissueRun = {t:0 for t in tissueEIDs}
	numRun = 0
	row = [0]*(3+numstates)
	row[0] = chm
	sql = "INSERT INTO epigenome VALUES (%s)" % (','.join('?' for n in xrange(3+numstates)))
	try:
		while True:
			chunk += 1
			row[1] = chunk
			for t in tissueEIDs:
				row[2] = t
				n = len(tissueEIDs[t])
				for s in xrange(numstates):
					del stateProbs[s][:]
				for e in tissueEIDs[t]:
					prob = list(float(p) for p in next(epiFile[e]).strip().split())
					for s,p in enumerate(prob):
						stateProbs[s].append(p)
				for s in xrange(numstates):
					stateProbs[s].sort()
					if n % 2:
						row[stateCol[s]] = stateProbs[s][n / 2]
					else:
						row[stateCol[s]] = (stateProbs[s][n / 2] + stateProbs[s][n / 2 - 1]) / 2.0
				run = max(stateCol, key=lambda c: row[c])
				if tissueRun[t] == run:
					numRun += 1
				else:
					tissueRun[t] = run
				cursor.execute(sql, row)
#				print "tissue #%d has %d epis, state %d = %g (%s)" % (t,n,numstates,row[stateCol[numstates-1]],str(stateProbs[numstates-1])) #TODO
			#for t
#			sys.exit(1) #TODO
		#while no StopIteration
	except StopIteration:
		pass
	
	# close all data files
	for fileobj in epiFile.itervalues():
		fileobj.close()
	epiFile = None
	
	cursor.execute("COMMIT TRANSACTION")
	print "  %d chunks (%d rows, %d runs)" % (chunk,chunk*len(tissueEIDs),numRun)
#for chm
