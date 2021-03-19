#!/usr/bin/env python

import argparse
import array
import bisect
import collections
import gc
import itertools
import resource
import struct
import sys
from zfile import zopen


class LDProfile:
	
	
	##################################################
	# class interrogation
	
	
	@classmethod
	def getVersionTuple(cls):
		# tuple = (major,minor,revision,dev,build,date)
		# dev must be in ('a','b','rc','release') for lexicographic comparison
		return (2,0,0,'a',3,'2014-03-07')
	#getVersionTuple()
	
	
	@classmethod
	def getVersionString(cls):
		v = list(cls.getVersionTuple())
		# tuple = (major,minor,revision,dev,build,date)
		# dev must be > 'rc' for releases for lexicographic comparison,
		# but we don't need to actually print 'release' in the version string
		v[3] = '' if v[3] > 'rc' else v[3]
		v[4] = v[4] or ''
		return "%d.%d.%d%s%s (%s)" % tuple(v)
	#getVersionString()
	
	
	##################################################
	# constructor
	
	
	def __init__(self):
		self._population = None
		self._chmFilter = None
		self._chmRefMap = None
		self._chmLabel = dict()
		self._chmRefPosition = dict()
		self._chmRefMarker = dict()
		self._chmRefOffset = dict()
		self._chmRefUprange = dict()
		self._chmRefDownrange = dict()
		self._chmBinaryDP = dict()
		self._chmBinaryRS = dict()
	#__init__()
	
	
	##################################################
	# public instance methods
	
	
	def getChromosomeFilter(self):
		return self._chmFilter
	#getChromosomeFilter()
	
	
	def setChromosomeFilter(self, chmFilter):
		self._chmFilter = set(c.strip().lower() for c in chmFilter) if chmFilter else None
	#setChromosomeFilter()
	
	
	def loadFromBinary(self, dataFile):
		# file header: magic bytes, format version number
		fmt = ">8sB"
		header,version = struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))
		if header != "RLAB_LDB":
			raise Exception("ERROR: missing LD binary format file header")
		
		if version == 1:
			# v1 file subheader: population label, number of chromosomes
			fmt = ">B"
			pop = dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0])
			if pop != (self._population or pop):
				raise Exception("ERROR: cannot mix data for populations %s and %s" % (self._population,pop))
			numChm = struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]
			
			# chromosome data chunks
			while numChm > 0:
				numChm -= 1
				
				# chunk header: chromosome label, number of reference points
				fmt = ">B"
				chm = dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0])
				fmt = ">I"
				numRef = struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]
				
				# skip filtered chromosomes
				if self._chmFilter and (chm.strip().lower() not in self._chmFilter):
					fmt = ">L"
					dataFile.seek(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0], 1)
					dataFile.seek(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0], 1)
					dataFile.seek(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0], 1)
					dataFile.seek(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0], 1)
					dataFile.seek(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0], 1)
					continue
				
				# initialize buffers
				refPosition = array.array('l')
				refMarker = None
				refOffset = array.array('l')
				refUprange = array.array('h')
				refDownrange = array.array('h')
				binaryDP = array.array('B')
				binaryRS = array.array('B')
				
				# reference point data
				fmt = ">L"
				refPosition.fromstring(dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]))
				if len(refPosition) != numRef:
					raise Exception("ERROR: unexpected reference point count on chromosome %s" % (chm,))
				refMarker = dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]).split("\0")
				refDownrange.fromstring(dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]))
				
				# regenerate parallel arrays
				refUprange.extend(0 for r in refDownrange)
				o = 0
				for r1,downrange in enumerate(refDownrange):
					for r2 in xrange(r1+1, r1+1+downrange):
						refUprange[r2] += 1
					refOffset.append(o)
					o += downrange
				
				# binary measurement data
				fmt = ">L"
				binaryDP.fromstring(dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]))
				if len(binaryDP) <= 0:
					binaryDP = None
				binaryRS.fromstring(dataFile.read(struct.unpack(fmt, dataFile.read(struct.calcsize(fmt)))[0]))
				if len(binaryRS) <= 0:
					binaryRS = None
				
				# audit and store the finished data structures
				self._finalizeLoad(pop, chm, refPosition, refMarker, refOffset, refUprange, refDownrange, binaryDP, binaryRS)
			#while numChm
		else:
			raise Exception("ERROR: unhandled file format version %d" % (version,))
		#if version
	#loadFromBinary()
	
	
	def dumpToBinary(self):
		# file header: magic bytes, format version number
		yield struct.pack(">8sB", "RLAB_LDB", 1)
		
		# v1 file subheader: population label, number of chromosomes
		chms = set(self._chmLabel)
		if self._chmFilter:
			chms &= self._chmFilter
		yield struct.pack(">B%dsB" % (len(self._population),),
			len(self._population), # under 2**8 (1 byte unsigned)
			self._population,
			len(chms) # under 2**8 (1 byte unsigned)
		)
		
		# chromosome data chunks
		for chm in sorted(chms):
			# chunk header: chromosome label, number of reference points
			label = self._chmLabel[chm]
			yield struct.pack(">B%dsI" % (len(label),),
				len(label), # under 2**8 (1 byte unsigned)
				label,
				len(self._chmRefPosition[chm])
			)
			
			# reference point data
			data = self._chmRefPosition[chm].tostring()
			yield struct.pack(">L%ds" % (len(data),), len(data), data)
			data = "\0".join(self._chmRefMarker[chm])
			yield struct.pack(">L%ds" % (len(data),), len(data), data)
			data = self._chmRefDownrange[chm].tostring()
			yield struct.pack(">L%ds" % (len(data),), len(data), data)
			
			# binary measurement data
			data = self._chmBinaryDP[chm].tostring() if (self._chmBinaryDP[chm] != None) else ""
			yield struct.pack(">L%ds" % (len(data),), len(data), data)
			data = self._chmBinaryRS[chm].tostring() if (self._chmBinaryRS[chm] != None) else ""
			yield struct.pack(">L%ds" % (len(data),), len(data), data)
		#foreach chm
	#dumpToBinary()
	
	
	def loadFromHapMap(self, chm, lineGenerator, loadDP=True, loadRS=True):
		if self._chmFilter and (chm.strip().lower() not in self._chmFilter):
			return
		filepop = None
		positionRef = dict()
		refPosition = array.array('l')
		markerRef = dict()
		refMarker = list()
		refOffset = array.array('l')
		refUprange = array.array('h')
		refDownrange = array.array('h')
		binaryDP = array.array('B') if loadDP else None
		binaryRS = array.array('B') if loadRS else None
		l = o = 0
		r1prev = r2prev = -1
		for line in lineGenerator:
			l += 1
			words = line.split(None,7)
			p1 = long(words[0])
			p2 = long(words[1])
			pop = words[2]
			m1 = words[3].strip().lower()
			m2 = words[4].strip().lower()
			dp = int(float(words[5]) / 0.005)
			rs = int(float(words[6]) / 0.005)
			
			# quick sanity checks
			if p1 >= p2:
				raise Exception("ERROR: columns 1 and 2 out of order on line %d: %d >= %d" % (l,p1,p2))
			elif (dp < 0) or (dp > 200):
				raise Exception("ERROR: column 6 invalid on line %d: %s" % (l,words[5]))
			elif (rs < 0) or (rs > 200):
				raise Exception("ERROR: column 7 invalid on line %d: %s" % (l,words[6]))
			elif filepop == None:
				if pop != (self._population or pop):
					raise Exception("ERROR: cannot mix data for populations %s and %s" % (self._population,pop))
				filepop = pop
			elif pop != filepop:
				raise Exception("ERROR: column 3 inconsistency on line %d: %s != %s" % (l,pop,filepop))
			
			# position/marker 1 sanity checks
			if p1 in positionRef:
				r1 = positionRef[p1]
				if m1 != refMarker[r1]:
					raise Exception("ERROR: wrong marker for position %d on line %d: %s != %s" % (p1,l,m1,refMarker.get(r1,'<None>')))
			elif refPosition and (p1 <= refPosition[-1]):
				raise Exception("ERROR: column 1 out of order on line %d: %d <= %d" % (l,p1,refPosition[-1]))
			elif m1 in markerRef:
				raise Exception("ERROR: wrong position for marker %s on line %d: %d != %d" % (m1,l,p1,refPosition[markerRef[m1]]))
			else:
				r1 = len(refPosition)
				positionRef[p1] = r1
				refPosition.append(p1)
				markerRef[m1] = r1
				refMarker.append(m1)
				refOffset.append(-1)
				refUprange.append(0)
				refDownrange.append(0)
			
			# position/marker 2 sanity checks
			if p2 in positionRef:
				r2 = positionRef[p2]
				if m2 != refMarker[r2]:
					raise Exception("ERROR: wrong marker for position %d on line %d: %s != %s" % (p2,l,m2,refMarker.get(r2,'<None>')))
			elif refPosition and (p2 <= refPosition[-1]):
				raise Exception("ERROR: column 2 out of order on line %d: %d <= %d" % (l,p2,refPosition[-1]))
			elif m2 in markerRef:
				raise Exception("ERROR: wrong position for marker %s on line %d: %d != %d" % (m2,l,p2,refPosition[markerRef[m2]]))
			else:
				r2 = len(refPosition)
				positionRef[p2] = r2
				refPosition.append(p2)
				markerRef[m2] = r2
				refMarker.append(m2)
				refOffset.append(-1)
				refUprange.append(0)
				refDownrange.append(0)
			
			# store the data values
			if r1 == r1prev:
				if r2 == r2prev + 1:
					# regular right-hand increment
					r2prev = r2
					refDownrange[r1] += 1
					assert (refDownrange[r1] > 1)
					refUprange[r2] += 1
					o += 1
					if loadDP:
						binaryDP.append(dp)
						assert (len(binaryDP) == o)
					if loadRS:
						binaryRS.append(rs)
						assert (len(binaryRS) == o)
				else:
					raise Exception("ERROR: columns 2/5 out of order on line %d: %d/%s != %d/%s" % (l,refPosition[r2],refMarker[r2],refPosition[r2prev+1],refMarker[r2prev+1]))
			elif r1 == r1prev + 1:
				if r2 == r1prev + 2:
					# regular left-hand transition
					r1prev = r1
					r2prev = r2
					assert (refOffset[r1] == -1)
					refOffset[r1] = o
					refDownrange[r1] += 1
					assert (refDownrange[r1] == 1)
					refUprange[r2] += 1
					o += 1
					if loadDP:
						binaryDP.append(dp)
						assert (len(binaryDP) == o)
					if loadRS:
						binaryRS.append(rs)
						assert (len(binaryRS) == o)
				else:
					raise Exception("ERROR: columns 2/5 out of order on line %d: %d/%s != %d/%s" % (l,refPosition[r2],refMarker[r2],refPosition[r1prev+2],refMarker[r1prev+2]))
			elif r1 == r1prev + 2:
				if r2 == r1prev + 3:
					# data gap transition
					r1prev = r1
					r2prev = r2
					assert (refOffset[r1-1] == -1)
					refOffset[r1-1] = o
					assert (refDownrange[r1-1] == 0)
					assert (refOffset[r1] == -1)
					refOffset[r1] = o
					refDownrange[r1] += 1
					assert (refDownrange[r1] == 1)
					refUprange[r2] += 1
					o += 1
					if loadDP:
						binaryDP.append(dp)
						assert (len(binaryDP) == o)
					if loadRS:
						binaryRS.append(rs)
						assert (len(binaryRS) == o)
				else:
					raise Exception("ERROR: columns 2/5 out of order on line %d: %d/%s != %d/%s" % (l,refPosition[r2],refMarker[r2],refPosition[r1prev+3],refMarker[r1prev+3]))
			else:
				raise Exception("ERROR: columns 1/4 out of order on line %d: %d/%s != %d/%s , %d/%s , %d/%s" % (l,refPosition[r1],refMarker[r1],refPosition[r1prev],refMarker[r1prev],refPosition[r1prev+1],refMarker[r1prev+1],refPosition[r1prev+2],refMarker[r1prev+2]))
		#for line in ldfile
		assert (refOffset[-1] == -1)
		refOffset[-1] = o
		
		# audit and store the finished data structures
		self._finalizeLoad(pop, chm, refPosition, refMarker, refOffset, refUprange, refDownrange, binaryDP, binaryRS)
	#loadFromHapMap()
	
	
	def dumpToHapMap(self, chm):
		if self._chmFilter and (chm.strip().lower() not in self._chmFilter):
			return
		pop = self._population
		refPosition = self._chmRefPosition[chm]
		refMarker = self._chmRefMarker[chm]
		refOffset = self._chmRefOffset[chm]
		refDownrange = self._chmRefDownrange[chm]
		binaryDP = self._chmBinaryDP[chm]
		binaryRS = self._chmBinaryRS[chm]
		for r1 in xrange(len(refPosition)):
			r2 = r1
			for o in xrange(refOffset[r1], refOffset[r1]+refDownrange[r1]):
				r2 += 1
				yield "%d %d %s %s %s %g %g\n" % (
					refPosition[r1], refPosition[r2],
					pop,
					refMarker[r1], refMarker[r2],
					binaryDP[o] * 0.005, binaryRS[o] * 0.005
				)
	#dumpToHapMap()
	
	
	def prune(self, minDP=0, minRS=0, reverse=False):
		minDP = min(max(0, int(minDP/0.005)), 200)
		minRS = min(max(0, int(minRS/0.005)), 200)
		
		for chm in self._chmLabel:
			# fetch original buffers
			oldRefPosition = self._chmRefPosition[chm]
			oldRefMarker = self._chmRefMarker[chm]
			oldRefOffset = self._chmRefOffset[chm]
			oldRefDownrange = self._chmRefDownrange[chm]
			oldBinaryDP = self._chmBinaryDP[chm]
			oldBinaryRS = self._chmBinaryRS[chm]
			refs = range(len(oldRefOffset))
			
			# initialize pruned buffers
			newRefOffset = array.array('l')
			newRefUprange = array.array('h')
			newRefDownrange = array.array('h')
			newBinaryDP = array.array('B')
			newBinaryRS = array.array('B')
			
			# scan and prune
			prevD = 0
			newO = 0
			for r in refs:
				oldO = oldRefOffset[r]
				oldD = oldRefDownrange[r]
				
				# set the new downrange to include the first below-threshold(s) value
				# (and don't let it decrease by more than 1, lest we create up/downstream gaps)
				if reverse:
					# decrease the downstream range as long as we're still below either threshold
					newD = oldD
					while (newD > 0) and ((oldBinaryRS[oldO+newD-1] < minRS) or (oldBinaryDP[oldO+newD-1] < minDP)):
						newD -= 1
					newD += 1
				else:
					# increase the downstream range as long as we're above both thresholds
					newD = 1
					while (newD < oldD) and (oldBinaryRS[oldO+newD-1] >= minRS) and (oldBinaryDP[oldO+newD-1] >= minDP):
						newD += 1
				#if reverse
				newD = min(max(0, prevD-1, newD), oldD)
				
				# copy over the un-pruned binary data
				newRefOffset.append(newO)
				newRefDownrange.append(newD)
				newBinaryDP.extend(oldBinaryDP[o] for o in xrange(oldO,oldO+newD))
				newBinaryRS.extend(oldBinaryRS[o] for o in xrange(oldO,oldO+newD))
				newO += newD
				prevD = newD
			#foreach ref
			
			# regenerate parallel arrays
			newRefUprange.extend(0 for r in refs)
			for r1,d in enumerate(newRefDownrange):
				for r2 in xrange(r1+1, r1+1+d):
					newRefUprange[r2] += 1
			
			# audit and store the pruned data structures
			self._finalizeLoad(self._population, chm, oldRefPosition, oldRefMarker, newRefOffset, newRefUprange, newRefDownrange, newBinaryDP, newBinaryRS)
		#foreach chm
	#prune()
	
	
	def getPopulation(self):
		return self._population
	#getPopulation()
	
	
	def getNumChromosomes(self):
		return len(self._chmRefPosition)
	#getNumChromosomes()
	
	
	def getChromosomes(self):
		return sorted(self._chmLabel.values())
	#getChromosome()
	
	
	def getNumPositions(self, chm):
		chm = str(chm).strip().lower()
		return len(self._chmRefPosition.get(chm,[]))
	#getNumPositions()
	
	
	def getNthPosition(self, chm, n):
		chm = str(chm).strip().lower()
		l = self._chmRefPosition.get(chm,[])
		if (n < -len(l)) or (n >= len(l)):
			return None
		return l[n]
	#getNthPosition()
	
	
	def enumeratePositions(self, chm, start=0):
		chm = str(chm).strip().lower()
		l = self._chmRefPosition.get(chm,[])
		if start < 0:
			return ( (r,l[r]) for r in xrange(len(l)+start, len(l)) )
		return ( (r,l[r]) for r in xrange(start, len(l)) )
	#enumeratePositions()
	
	
	def getNumMarkers(self, chm):
		chm = str(chm).strip().lower()
		return len(self._chmRefMarker.get(chm,[]))
	#getNumMarkers()
	
	
	def getNthMarker(self, chm, n):
		chm = str(chm).strip().lower()
		l = self._chmRefMarker.get(chm,[])
		if (n < -len(l)) or (n >= len(l)):
			return None
		return l[n]
	#getNthMarker()
	
	
	def enumerateMarkers(self, chm, start=0):
		chm = str(chm).strip().lower()
		l = self._chmRefMarker.get(chm,[])
		if start < 0:
			return ( (r,l[r]) for r in xrange(len(l)+start, len(l)) )
		return ( (r,l[r]) for r in xrange(start, len(l)) )
	#enumerateMarkers()
	
	
	def enumerateMarkerPositions(self, chm, start=0):
		chm = str(chm).strip().lower()
		lm = self._chmRefMarker.get(chm,[])
		lp = self._chmRefPosition.get(chm,[])
		if start < 0:
			return ( (r,lm[r],lp[r]) for r in xrange(len(lm)+start, len(lm)) )
		return ( (r,lm[r],lp[r]) for r in xrange(start, len(lm)) )
	#enumerateMarkerPositions()
	
	
	def getNumDPValues(self, chm):
		chm = str(chm).strip().lower()
		return len(self._chmBinaryDP.get(chm,""))
	#getNumDPValues()
	
	
	def getDP(self, chm, pos1, pos2):
		chm = str(chm).strip().lower()
		return self._getInterpolatedValue(self._chmBinaryDP[chm], chm, pos1, pos2)
	#getDP()
	
	
	def auditDP(self, chm, skipUp=True, skipDown=True):
		chm = str(chm).strip().lower()
		return self._audit(chm, self._chmBinaryDP[chm], skipUp, skipDown)
	#auditDP()
	
	
	def getNumRSValues(self, chm):
		chm = str(chm).strip().lower()
		return len(self._chmBinaryRS.get(chm,""))
	#getNumRSValues()
	
	
	def getRS(self, chm, pos1, pos2):
		chm = str(chm).strip().lower()
		return self._getInterpolatedValue(self._chmBinaryRS[chm], chm, pos1, pos2)
	#getRS()
	
	
	def auditRS(self, chm, skipUp=True, skipDown=True):
		chm = str(chm).strip().lower()
		return self._audit(chm, self._chmBinaryRS[chm], skipUp, skipDown)
	#auditRS()
	
	
	##################################################
	# private instance methods
	
	
	def _finalizeLoad(self, pop, chm, refPosition, refMarker, refOffset, refUprange, refDownrange, binaryDP, binaryRS):
		# validate population consistency and length
		assert (pop and (len(pop) < 2**8))
		assert (pop == (self._population or pop))
		
		# validate chromosome count and label length
		label = str(chm)
		chm = label.strip().lower()
		assert ((len(self._chmLabel) + (0 if (chm in self._chmLabel) else 1)) < 2**8)
		assert (len(label) < 2**8)
		
		# make sure all reference property arrays are the same length
		assert (len(refPosition) == len(refMarker) == len(refOffset) == len(refDownrange))
		
		# validate reference point count
		numRef = len(refPosition)
		assert (numRef < 2**32)
		
		# validate reference position values
		assert ((min(refPosition) >= 0) and (max(refPosition) < 2**63))
		
		# make sure the reference positions are in order with no duplicates
		assert (all( (refPosition[r-1] < refPosition[r]) for r in xrange(1,numRef) ))
		
		# validate total marker label length
		assert (sum(len(m) for m in refMarker) < 2**32)
		
		# make sure the marker list has no duplicates
		assert (len(set(refMarker)) == numRef)
		
		# validate up/downstream range values
		assert ((min(refDownrange) >= 0) and (max(refDownrange) < 2**15))
		assert ((min(refUprange) >= 0) and (max(refUprange) < 2**15))
		
		# make sure the downrange tallies agree with the uprange tallies and binary offsets and sizes
		o = d1min = u1max = 0
		for r1 in xrange(numRef):
			u1 = refUprange[r1]
			d1 = refDownrange[r1]
			assert (o == refOffset[r1])
			assert (d1 >= d1min)
			assert (u1 <= u1max)
			r2 = r1 + d1
			if r2 < len(refUprange):
				assert (refUprange[r2] >= d1)
			r2 += 1
			if r2 < len(refUprange):
				assert (refUprange[r2] <= d1)
			d1min = d1 - 1
			u1max = u1 + 1
			o += d1
		assert (o == sum(refUprange))
		assert ((binaryDP == None) or (o == len(binaryDP)))
		assert ((binaryRS == None) or (o == len(binaryRS)))
		
		# validate binary data size
		assert (o < 2**64)
		
		# store everything on the instance
		self._population = self._population or pop
		self._chmLabel[chm] = label
		self._chmRefPosition[chm] = refPosition
		self._chmRefMarker[chm] = refMarker
		self._chmRefOffset[chm] = refOffset
		self._chmRefUprange[chm] = refUprange
		self._chmRefDownrange[chm] = refDownrange
		self._chmBinaryDP[chm] = binaryDP
		self._chmBinaryRS[chm] = binaryRS
	#_finalizeLoad()
	
	
	def _getRefWeights(self, chm, pos, skip1=None, skip2=None):
		chm = str(chm).strip().lower()
		refPosition = self._chmRefPosition.get(chm,[])
		
		rR = bisect.bisect_left(refPosition, pos)
		rL = rR - 1
		while (skip1 == rL) or (skip2 == rL):
			rL -= 1
		while (skip1 == rR) or (skip2 == rR):
			rR += 1
		
		if rR >= len(refPosition):
			if rL < 0:
				return tuple()
			return ( (rL, 1.0), )
		elif (rL < 0) or (pos == refPosition[rR]):
			return ( (rR, 1.0), )
		wR = float(pos - refPosition[rL]) / (refPosition[rR] - refPosition[rL])
		return ( (rL, 1.0-wR), (rR, wR) )
	#_getRefWeights()
	
	
	def _getInterpolatedValue(self, binary, chm, pos1, pos2, skip1=None, skip2=None):
	#	print "LD%s(%d, %d) =" % (chm,pos1,pos2)
		chm = str(chm).strip().lower()
		refOffset = self._chmRefOffset[chm]
		refDownrange = self._chmRefDownrange[chm]
		weights1 = self._getRefWeights(chm, pos1, skip1, skip2)
		weights2 = self._getRefWeights(chm, pos2, skip1, skip2)
		w = v = 0.0
		for r1,w1 in weights1:
			for r2,w2 in weights2:
	#			print "  %1.2f * LD%s(%d, %d) * %1.2f *" % (w1, chm, r1, r2, w2),
				w += w1 * w2
				if r1 == r2:
	#				print "%1.2f" % (1.0,)
					v += w1 * w2 # * 1.0
				else:
					if r1 > r2:
						r1,r2 = r2,r1
					if r2 <= r1 + refDownrange[r1]:
	#					print "%1.2f" % (binary[refOffset[r1] + r2 - r1 - 1] * 0.005,)
						v += w1 * w2 * binary[refOffset[r1] + r2 - r1 - 1] * 0.005
					#else += 0.0
		return ( (v / w), w )
	#_getInterpolatedValue()
	
	
	def _audit(self, chm, binary, skipUp=True, skipDown=True):
		chm = str(chm).strip().lower()
		refPosition = self._chmRefPosition[chm]
		refUprange = self._chmRefUprange[chm]
		refDownrange = self._chmRefDownrange[chm]
		refOffset = self._chmRefOffset[chm]
		if True: #TODO: test all points, even edge cases?
			for r1,p1 in enumerate(refPosition):
				e = list()
				for r2 in xrange(r1+1, r1+1+refDownrange[r1]):
					vAct = binary[refOffset[r1] + r2 - r1 - 1] * 0.005
					vEst,_ = self._getInterpolatedValue(binary, chm, p1, refPosition[r2], (r1 if skipUp else None), (r2 if skipDown else None))
					e.append(vEst - vAct)
				yield e
		elif skipUp and skipDown:
			for r1 in xrange(1,len(refPosition)-1):
				p1 = refPosition[r1]
				e = list()
				for r2 in xrange(r1+3, r1+refDownrange[r1]):
					vAct = binary[refOffset[r1] + r2 - r1 - 1] * 0.005
					vEst,_ = self._getInterpolatedValue(binary, chm, p1, refPosition[r2], r1, r2)
					e.append(vEst - vAct)
				yield e
		elif skipUp:
			for r1 in xrange(1,len(refPosition)-1):
				p1 = refPosition[r1]
				e = list()
				for r2 in xrange(r1+2, r1+1+refDownrange[r1]):
					vAct = binary[refOffset[r1] + r2 - r1 - 1] * 0.005
					vEst,_ = self._getInterpolatedValue(binary, chm, p1, refPosition[r2], r1, None)
					e.append(vEst - vAct)
				yield e
		elif skipDown:
			for r1 in xrange(0,len(refPosition)):
				p1 = refPosition[r1]
				e = list()
				for r2 in xrange(r1+2, r1+1+refDownrange[r1]):
					vAct = binary[refOffset[r1] + r2 - r1 - 1] * 0.005
					vEst,_ = self._getInterpolatedValue(binary, chm, p1, refPosition[r2], None, r2)
					e.append(vEst - vAct)
				yield e
	#_audit()
	
	
#class LDProfile


if __name__ == "__main__":
	# zcat ld_chr22_CEU.txt.gz | awk '{print $1,$2,$3,$4,$5,0.005*int($6/0.005),0.005*int($7/0.005)}' > ceu.txt
	# zcat ld_chr2_YRI.txt.gz  | awk '{print $1,$2,$3,$4,$5,0.005*int($6/0.005),0.005*int($7/0.005)}' > yri.txt
	
	version = "LDSpline version %s" % (LDProfile.getVersionString(),)
	
	# define arguments
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=version,
	)
	parser.add_argument('--version', action='version', version=version)
	#TODO precision option? 1 byte (0.005), 2 byte (0.00002), 3 byte (1e-7), 4 byte (...)
	parser.add_argument('-c', '--chromosomes', type=str, metavar='chr', nargs='+', action='append', default=None,
			help="reduce memory usage by only loading LD data for the specified chromosomes",
	)
	parser.add_argument('-l', '--liftover-input', type=str, metavar=('oldfile','newfile','unmapfile'), nargs=3, action='store', default=None,
			help="use the given liftOver results to update the LD profile's marker map"
	)
	parser.add_argument('-b', '--binary-input', type=str, metavar='file', action='append', default=None,
			help="read LD profile data from the given file in LDSpline binary format"
	)
	parser.add_argument('-t', '--text-input', type=str, metavar=('chr','file'), nargs=2, action='append', default=None,
			help="read LD profile data for the given chromosome from the given file in HapMap text format"
	)
	parser.add_argument('--prune-dp', type=float, metavar='minDP', action='store', default=None,
			help="reduce the size of the LD profile by setting a minimum d' value and dropping data beyond that threshold",
	)
	parser.add_argument('--prune-rs', type=float, metavar='minRS', action='store', default=None,
			help="reduce the size of the LD profile by setting a minimum r^2 value and dropping data beyond that threshold",
	)
	parser.add_argument('--prune-more', action='store_true', default=False,
			help="prune more data by scanning outward instead of inward (default: off)",
	)
	parser.add_argument('--audit', action='store_true', default=None,
			help="audit the interpolation algorithm by comparing all data values to their interpolated estimates",
	)
	parser.add_argument('-M', '--marker-output', type=str, metavar='file', action='store', default=None,
			help="export all LD profile markers to the given file in 3-column .MAP format"
	)
	parser.add_argument('-L', '--liftover-output', type=str, metavar='file', action='store', default=None,
			help="export all LD profile markers to the given file in .BIM format, suitable for liftOver"
	)
	parser.add_argument('-B', '--binary-output', type=str, metavar='file', action='store', default=None,
			help="write LD profile data to the given file in LDSpline binary format"
	)
	parser.add_argument('-T', '--text-output', type=str, metavar=('chr','file'), nargs=2, action='append', default=None,
			help="write LD profile data for the given chromosome to the given file in HapMap text format"
	)
	
	# if no arguments, print usage and exit
	if len(sys.argv) < 2:
		print version
		print
		parser.print_usage()
		print
		print "Use -h for details."
		sys.exit(2)
	
	# initialize and parse arguments
	ldp = LDProfile()
	args = parser.parse_args()
	empty = tuple()
	
	# set chromosome filter
	chms = sorted(set(itertools.chain(*(args.chromosomes or empty))))
	if chms:
		print "setting chromosome filter: %s" % (', '.join(chms),)
		ldp.setChromosomeFilter(chms)
	
	# read liftOver results?
	if args.liftover_input:
		loChrOld = collections.defaultdict(list)
		loChrNew = collections.defaultdict(list)
		loChrDrop = collections.defaultdict(list)
		print "reading liftOver input from '%s' ..." % (args.liftover_input[0],)
		n = 0
		with open(args.liftover_input[0],'rU') as datafile:
			for line in datafile:
				n += 1
				c,p,_ = line.split(None,2)
				c = c.lower()
				if c.startswith('chr'):
					c = c[3:]
				loChrOld[c].append(long(p))
		print "... OK: read %d markers on %d chromosomes" % (n,len(loChrOld))
		print "reading liftOver updates from '%s' ..." % (args.liftover_input[1],)
		n = 0
		with open(args.liftover_input[1],'rU') as datafile:
			for line in datafile:
				n += 1
				c,p,_ = line.split(None,2)
				c = c.lower()
				if c.startswith('chr'):
					c = c[3:]
				loChrNew[c].append(long(p))
		print "... OK: read %d markers on %d chromosomes" % (n,len(loChrNew))
		print "reading liftOver drops from '%s' ..." % (args.liftover_input[2],)
		n = 0
		with open(args.liftover_input[2],'rU') as datafile:
			for line in datafile:
				if not line.startswith('#'):
					n += 1
					c,p,_ = line.split(None,2)
					c = c.lower()
					if c.startswith('chr'):
						c = c[3:]
					loChrDrop[c].append(long(p))
		print "... OK: read %d markers on %d chromosomes" % (n,len(loChrDrop))
		print "analyzing liftOver results ..."
		loChrMap = collections.defaultdict(dict)
		for c,lOld in loChrOld.iteritems():
			xOld = len(lOld)
			lNew = loChrNew.get(c,empty)
			xNew = len(lNew)
			lDrop = loChrDrop.get(c,empty)
			xDrop = len(lDrop)
			nOld = nNew = nDrop = 0
			revMap = dict()
			while nOld < xOld:
				pOld = lOld[nOld]
				if (nDrop < xDrop) and (lOld[nOld] == lDrop[nDrop]):
					pNew = None
					nDrop += 1
				elif (nNew < xNew):
					pNew = lNew[nNew]
					nNew += 1
				else:
					print "ERROR: liftOver record count mismatch at chr%s:%d" % (c,pOld)
					exit(1)
				
				if pNew and (revMap.get(pNew,pOld) != pOld):
		#			print "ERROR: liftOver result conflation for chr%s:%d <- %d , %d" % (c,pNew,revMap[pNew],pOld)
		#			exit(1)
					pNew = None
				if loChrMap[c].get(pOld,pNew) != pNew:
					print "ERROR: liftOver result mismatch for chr%s:%d -> %s , %s" % (c,pOld,loChrMap[c][pOld],pNew)
					exit(1)
				loChrMap[c][pOld] = pNew
				if pNew:
					revMap[pNew] = pOld
				
				nOld += 1
			#while nOld < xOld
			if nNew < xNew:
				print "WARNING: excess liftOver updates on chr%s" % (c,)
			if nDrop < xDrop:
				print "WARNING: excess liftOver drops on chr%s" % (c,)
		#foreach loChrOld
		loChrOld = loChrNew = loChrDrop = None
		print "... OK"
	#if liftover_input
	
	# read binary data?
	for filename in (args.binary_input or empty):
		print "reading LD profile data from '%s' in binary format ..." % (filename,)
		with open(filename,'rb') as datafile:
			ldp.loadFromBinary(datafile)
		print "... OK: profile covers %d chromosome(s), %d referenece points, %d measurements" % (
			ldp.getNumChromosomes(),
			sum(ldp.getNumMarkers(c) for c in ldp.getChromosomes()),
			sum(ldp.getNumDPValues(c) for c in ldp.getChromosomes())
		)
	#foreach binary_input
	
	# read text data?
	for chm,filename in (args.text_input or empty):
		print "reading chr %s LD profile data from '%s' in text format ..." % (chm,filename)
		with zopen(filename) as datafile:
			ldp.loadFromHapMap(chm, datafile)
		print "... OK: profile covers %d chromosome(s), %d referenece points, %d measurements" % (
			ldp.getNumChromosomes(),
			sum(ldp.getNumMarkers(c) for c in ldp.getChromosomes()),
			sum(ldp.getNumDPValues(c) for c in ldp.getChromosomes())
		)
	#foreach text_input
	
	
	# prune?
	if args.prune_dp or args.prune_rs:
		print "pruning LD profile to d' >= %g and r^2 >= %g ..." % ((args.prune_dp or 0.0), (args.prune_rs or 0.0))
		d0 = sum(ldp.getNumRSValues(c) for c in ldp.getChromosomes())
		ldp.prune((args.prune_dp or 0.0), (args.prune_rs or 0.0), not args.prune_more)
		d1 = sum(ldp.getNumRSValues(c) for c in ldp.getChromosomes())
		print "... OK: dropped %d data points (%d%% reduction)" % (d0-d1,100*(d0-d1)/d0)
	#if prune
	
	# audit?
	if args.audit:
		for metric,method in (("d'",ldp.auditDP),("r^2",ldp.auditRS)):
			print "auditing %s interpolation error ..." % (metric,)
			tn = te = te2 = 0.0
			for case,up,dn in (('upstream',1,0),('downstream',0,1),('both',1,1)):
				n = e = e2 = 0.0
				for c in ldp.getChromosomes():
					for referrors in method(c, up, dn):
						n += len(referrors)
						for error in referrors:
							e += abs(error)
							e2 += error**2
				print "  %s: MAE=%1.3f, RMSE=%1.3f (over %d points)" % (case,e/n,(e2/n)**0.5,n)
				tn += n
				te += e
				te2 += e2
			print "... OK: total MAE=%1.3f, RMSE=%1.3f (over %d points)" % (te/tn,(te2/tn)**0.5,tn)
	#if audit
	
	
	# write marker map?
	if args.marker_output:
		print "writing markers to '%s' in 3-column .MAP format ..." % (args.marker_output,)
		l = 0
		with open(args.marker_output,'wb') as datafile:
			for c in ldp.getChromosomes():
				for n,m,p in ldp.enumerateMarkerPositions(c):
					l += 1
					datafile.write("%s\t%s\t%d\n" % (c,m,p))
		print "... OK: wrote %d markers" % (l,)
	#if marker_output
	
	# write markers for liftover?
	if args.liftover_output:
		print "writing markers to '%s' in .BIM format (for liftOver) ..." % (args.liftover_output,)
		l = 0
		with open(args.liftover_output,'wb') as datafile:
			for c in ldp.getChromosomes():
				for n,p in ldp.enumeratePositions(c):
					l += 1
					datafile.write("chr%s\t%d\t%d\n" % (c,p,p+1))
		print "... OK: wrote %d markers" % (l,)
	#if liftover_output
	
	# write binary data?
	if args.binary_output:
		print "writing LD profile data to '%s' in binary format ..." % (args.binary_output,)
		b = 0
		with open(args.binary_output,'wb') as datafile:
			for chunk in ldp.dumpToBinary():
				datafile.write(chunk)
				b += len(chunk)
		print "... OK: wrote %d bytes" % (b,)
	#if binary_output
	
	# write text data?
	for chm,filename in (args.text_output or empty):
		print "writing chr %s LD profile data to '%s' in text format ..." % (chm,filename)
		l = 0
		with open(filename,'wb') as datafile:
			for line in ldp.dumpToHapMap(chm):
				datafile.write(line)
				l += 1
		print "... OK: wrote %d lines" % (l,)
	#foreach text_output
	
#if __main__
