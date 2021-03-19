#!/usr/bin/env python

import struct
import zlib


def TYPE_UNKNOWN = 0
def TYPE_ASCII = 1
def TYPE_ZLIB = 2
def TYPE_GZIP = 3
def TYPE_BZIP = 4


class zfile(object):
	
	
	def __init__(self, filename, mode=None, compresstype=None, compresslevel=None, fileobj=None, mtime=None, newlines="\n", bufsize=512*1024):
		_self = self.__dict__
		
		_self['_file'] = fileobj or open(filename, mode or 'rb')
		_self['_name'] = filename or _self['_file'].name
		_self['_mode'] = mode.lower()
		_self['_type'] = compresstype
		_self['_level'] = compresslevel
		_self['_mtime'] = mtime
		_self['_newlines'] = newlines
		_self['_bufsize'] = max(1,bufsize)
		
		_self['_readonly'] = set(['closed','encoding','errors','mode','name','newlines'])
		_self['_closed'] = False
		_self['_encoding'] = None
		_self['_errors'] = None
		
		_self['_readbuffer'] = ''
		_self['_textbuffer'] = ''
		_self['_linebuffer'] = []
		_self['_writebuffer'] = ''
		
		_self['_zlib_dc'] = None
		
		_self['_gzip_ftext'] = None
		_self['_gzip_fhcrc'] = None
		_self['_gzip_fextra'] = None
		_self['_gzip_fname'] = None
		_self['_gzip_fcomment'] = None
		_self['_gzip_mtime'] = None
		_self['_gzip_xflslow'] = None
		_self['_gzip_xflfast'] = None
		_self['_gzip_os'] = None
		_self['_gzip_crc32'] = None
		_self['_gzip_isize'] = None
	#__init__()
	
	
	def __del__(self):
		self.close()
	#__del__()
	
	
	def close(self):
		_self = self.__dict__
		if not _self['_closed']:
			_self['_closed'] = True
			_self['_file'].close()
	#close()
	
	
	def flush(self):
		pass #TODO
	#flush()
	
	
	def fileno(self):
		return self.__dict__['_file'].fileno()
	#fileno()
	
	
	def isatty(self):
		return self.__dict__['_file'].isatty()
	#isatty()
	
	
	def next(self):
		return self.__next__()
	#next()
	
	
	def _rawreaduntil(self, size):
		_self = self.__dict__
		while len(_self['_readbuffer']) < size:
			data = _self['_file'].read(_self['_bufsize'])
			if not data:
				raise StopIteration
			_self['_readbuffer'] += data
	#_rawreaduntil()
	
	
	def read(self, size=None):
		_self = self.__dict__
		
		# if we don't know the file type yet, peek at the first chunk and try to guess
		if not _self['_type']:
			_self['_type'] = TYPE_ASCII
			try:
				self._rawreaduntil(4)
				b0 = ord(_self['_readbuffer'][0])
				b1 = ord(_self['_readbuffer'][1])
				b2 = ord(_self['_readbuffer'][2])
				b3 = ord(_self['_readbuffer'][3])
				if (b0 == 0x1f) and (b1 == 0x8b) and (b2 == 0x08) and ((b3 & 0xe0) == 0):
					_self['_type'] = TYPE_GZIP
				elif (b0 & 0x0F == 8) and (((b0 & 0xF0) >> 4) <= 7) and (((b0*256 + b1) % 31) == 0):
					_self['_type'] = TYPE_ZLIB
			except StopIteration:
				pass
		#if type unknown
		
		# validate the file header
		if _self['_type'] == TYPE_GZIP:
			try:
				dpos = 10
				self._rawreaduntil(dpos)
				flg = ord(_self['_readbuffer'][3])
				ftext = (flg & 0x01) != 0
				fhcrc = (flg & 0x02) != 0
				fextra = (flg & 0x04) != 0
				fname = (flg & 0x08) != 0
				fcomment = (flg & 0x10) != 0
				_self['_gzip_mtime'] = struct.unpack('<i', _self['_readbuffer'][4:8])
				xfl = ord(_self['_readbuffer'][8])
				_self['_gzip_xflslow'] = (xfl & 0x02) != 0
				_self['_gzip_xflfast'] = (xfl & 0x04) != 0
				_self['_gzip_os'] = ord(_self['_readbuffer'][9])
				
				if fextra:
					self._rawreaduntil(dpos+2)
					xlen = struct.unpack('<w', _self['_readbuffer'][dpos:dpos+2])
					hpos,dpos = dpos+2,dpos+2+xlen
					self._rawreaduntil(dpos)
					_self['_gzip_fextra'] = _self['_readbuffer'][hpos:dpos]
				
				if fname:
					hpos,dpos = dpos,_self['_readbuffer'].find('\0',dpos)
					while dpos < 0:
						self._rawreaduntil(len(_self['_readbuffer']) + 1)
						dpos = _self['_readbuffer'].find('\0',hpos)
					_self['_gzip_fname'] = _self['_readbuffer'][hpos:dpos]
					dpos += 1
				
				if fcomment:
					hpos,dpos = dpos,_self['_readbuffer'].find('\0',dpos)
					while dpos < 0:
						self._rawreaduntil(len(_self['_readbuffer']) + 1)
						dpos = _self['_readbuffer'].find('\0',hpos)
					_self['_gzip_fcomment'] = _self['_readbuffer'][hpos:dpos]
					dpos += 1
				
				if fhcrc:
					self._rawreaduntil(dpos+2)
					crc16 = struct.unpack('<w', _self['_readbuffer'][dpos:dpos+2])
					dpos += 2
				
				_self['_gzip_crc32'] = None
				_self['_gzip_isize'] = None
				
				# strip off the gzip header and remap read() to use zlib decompression
				_self['_readbuffer'] = _self['_readbuffer'][dpos:]
				self.read = self._read_zlib
				return self._read_zlib(size)
				
			except StopIteration:
				_self['_type'] = TYPE_ASCII
			
		elif _self['_type'] == TYPE_ZLIB:
			self.read = self._read_zlib
			return self._read_zlib(size)
			
		else:
			self.read = self._read_ascii
			return self._read_ascii(size)
	#read()
	
	
	def _read_ascii(self, size=None):
		_self = self.__dict__
		# if split lines remain, re-merge them with unsplit buffer
		if _self['_linebuffer']:
			_self['_linebuffer'].append(_self['_readbuffer'])
			_self['_readbuffer'] = _self['_newlines'].join(_self['_linebuffer'])
			_self['_linebuffer'] = []
		# if no size was given, read and return everything
		if size <= 0:
			data = _self['_readbuffer'] + _self['_file'].read()
			_self['_readbuffer'] = ''
			return data
		# read until we have enough in the buffer
		try:
			self._rawreaduntil(size)
		except StopIteration:
			pass
		# return as much as requested and leave the rest in the buffer
		data = _self['_readbuffer'][0:size]
		_self['_readbuffer'] = _self['_readbuffer'][len(data):]
		return data
	#_read_ascii()
	
	
	def _read_zlib(self, size=None):
		_self = self.__dict__
		# if split lines remain, re-merge them with unsplit buffer
		if _self['_linebuffer']:
			_self['_linebuffer'].append(_self['_textbuffer'])
			_self['_textbuffer'] = _self['_newlines'].join(_self['_linebuffer'])
			_self['_linebuffer'] = []
		# create a decompressor if we don't have one yet
		if not _self['_zlib_dc']:
				_self['_zlib_dc'] = zlib.decompressobj(zlib.MAX_WBITS) # | 32)
		# if no size was given, read, decompress and return everything
		if size <= 0:
			data = _self['_textbuffer'] + _self['_zlib_dc'].decompress(_self['_readbuffer'] + _self['_file'].read()) + _self['_zlib_dc'].flush()
			_self['_readbuffer'] = ''
			_self['_textbuffer'] = ''
			_self['_zlib_dc'] = None
			return data
		# read and decompress until we have enough in the buffer
		while len(_self['_textbuffer']) < size:
			data = _self['_file'].read(_self['_bufsize'])
			if (not data) and (not _self['_readbuffer']):
				_self['_textbuffer'] += _self['_zlib_dc'].flush()
				_self['_zlib_dc'] = None
				break
			_self['_textbuffer'] += _self['_zlib_dc'].decompress(_self['_readbuffer'] + data)
			_self['_readbuffer'] = ''
		# return as much as requested and leave the rest in the buffer
		data = _self['_textbuffer'][0:size]
		_self['_textbuffer'] = _self['_textbuffer'][len(data):]
		return data
	#_read_zlib()
	
	
	def readline(self, size=None):
		pass #TODO
	#readline()
	
	
	def readlines(self, sizehint=None):
		pass #TODO
	#readlines()
	
	
	def seek(self, offset, whence=None):
		pass #TODO
	#seek()
	
	
	def tell(self):
		pass #TODO
	#tell()
	
	
	def truncate(self, size=None):
		pass #TODO
	#truncate()
	
	
	def write(self, data):
		pass #TODO
	#write()
	
	
	def writelines(self, sequence):
		pass #TODO
	#writelines()
	
	
	def __getattr__(self, name):
		_self = self.__dict__
		if name.startswith('_'):
			raise AttributeError("attribute '%s' of 'zfile' objects is not public" % name)
		if name in _self['_readonly']:
			return _self['_'+name]
		if name not in _self:
			raise AttributeError("'zfile' object has no attribute '%s'" % name)
		return _self[name]
	#__getattr__()
	
	
	def __setattr__(self, name, value):
		_self = self.__dict__
		if name.startswith('_'):
			raise AttributeError("attribute '%s' of 'zfile' objects is not public" % name)
		if name in _self['_readonly']:
			raise AttributeError("attribute '%s' of 'zfile' objects is not writable" % name)
		_self[name] = value
	#__setattr__()
	
	
	def __delattr__(self, name):
		_self = self.__dict__
		if name.startswith('_'):
			raise AttributeError("attribute '%s' of 'zfile' objects is not public" % name)
		if name in _self['_readonly']:
			raise AttributeError("attribute '%s' of 'zfile' objects is not writable" % name)
		if name not in _self:
			raise AttributeError("'zfile' object has no attribute '%s'" % name)
		del _self[name]
	#__delattr__()
	
	
	def __enter__(self):
		return self
	#__enter__()
	
	
	def __exit__(self, excType, excVal, excTrace):
		self.close()
	#__exit__()
	
	
	def __iter__(self):
		return self
	#__iter__()
	
	
	def __next__(self):
		#TODO
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
			return self.__next()
		# reverse the remaining lines into a stack and pop one to return
		self._lines.reverse()
		return self._lines.pop()
	#__next__()
	
#zfile


if __name__ == '__main__':
	import sys
	for n in xrange(1, len(sys.argv)):
		with zfile(sys.argv[n]) as f:
			for line in f:
				print line
