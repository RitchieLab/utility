#!/usr/bin/env python

import zlib


class zopen(object):
	
	def __init__(self, fileName, splitChar="\n", chunkSize=16*1024):
		self._filePtr = open(fileName,'rb')
		self._splitChar = splitChar
		self._chunkSize = chunkSize
		self._dc = zlib.decompressobj(zlib.MAX_WBITS+32) # autodetect gzip or zlib header
		self._text = ""
		self._lines = list()
	#__init__()
	
	
	def __del__(self):
		if self._filePtr:
			self._filePtr.close()
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
			data = self._filePtr.read(self._chunkSize)
			if data:
				self._text += self._dc.decompress(data)
				del data
			else:
				self._text += self._dc.flush()
				del self._dc
		# if there's no text left, we're done
		if not self._text:
			raise StopIteration
		# split the text into lines
		self._lines = self._text.split(self._splitChar)
		self._text = ""
		# if there's more than one line, store the last to combine with the next chunk
		if len(self._lines) > 1:
			self._text = self._lines.pop()
		# reverse the remaining lines into a stack and pop one to return
		self._lines.reverse()
		return self._lines.pop()
	#__next__()
	
	
	def next(self):
		return self.__next__()
	#next()
	
#zopen


if __name__ == '__main__':
	import sys
	for n in xrange(1, len(sys.argv)):
		with zopen(sys.argv[n]) as zfile:
			for line in zfile:
				print line
