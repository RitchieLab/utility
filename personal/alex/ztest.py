#!/usr/bin/env python

import gzip
import struct
import subprocess
import sys
import time
import zlib


if len(sys.argv) <= 2:
	sys.stderr.write("usage: %s <compress|decompress|lines> <pipe|zlib|gzip> [chunksize] [filename]\n")
	sys.exit(2)
chunkSize = int(sys.argv[3]) if len(sys.argv) > 3 else 0
chunkSize = chunkSize or 512*1024
fileName = sys.argv[4] if len(sys.argv) > 4 else None

if sys.argv[1] == 'compress':
	if sys.argv[2] == 'pipe':
		sys.stderr.write("baseline (os pipe) compression test...\n")
		t0 = time.time()
		subprocess.call(["gzip","-6","-n","-c"])
		t1 = time.time()
		sys.stderr.write("...completed in %f seconds\n" % (t1-t0))
		
	elif sys.argv[2] == 'zlib':
		sys.stderr.write("zlib module compression test (chunksize=%d)...\n" % chunkSize)
		t0 = time.time()
		# write RFC1952 gzip header
		sys.stdout.write("\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff") # gzip magic numbers, mode=deflate, flags=0, timestamp=0, xflags=0, os=?
		cobj = zlib.compressobj(6)
		crc32 = 0
		size = 0
		last = None
		data = sys.stdin.read(chunkSize)
		if data:
			# skip 2-byte RFC1950 zlib header (compression format, flags)
			crc32 = zlib.crc32(data)
			size = len(data)
			sys.stdout.write(cobj.compress(data)[2:])
			data = sys.stdin.read(chunkSize)
			while data:
				crc32 = zlib.crc32(data, crc32)
				size += len(data)
				sys.stdout.write(cobj.compress(data))
				data = sys.stdin.read(chunkSize)
		# skip 2-byte RFC1950 zlib footer (ADLTER32 checksum)
		sys.stdout.write(cobj.flush()[:-4])
		sys.stdout.write(struct.pack('I',crc32 & 0xffffffff)) # crc32 of uncompressed data
		sys.stdout.write(struct.pack('I',size & 0xffffffff)) # size of uncompressed data modulo 2^32
		t1 = time.time()
		sys.stderr.write("...completed in %f seconds\n" % (t1-t0))
		
	elif sys.argv[2] == 'gzip':
		sys.stderr.write("gzip module compression test (chunksize=%d)...\n" % chunkSize)
		t0 = time.time()
		with gzip.GzipFile('', 'wb', compresslevel=6, fileobj=sys.stdout, mtime=0) as f:
			data = sys.stdin.read(chunkSize)
			while data:
				f.write(data)
				data = sys.stdin.read(chunkSize)
		t1 = time.time()
		sys.stderr.write("...completed in %f seconds\n" % (t1-t0))
		
elif sys.argv[1] == 'decompress':
	if sys.argv[2] == 'pipe':
		sys.stderr.write("baseline (os pipe) decompression test...\n")
		t0 = time.time()
		subprocess.call(["gunzip","-c"])
		t1 = time.time()
		sys.stderr.write("...completed in %f seconds\n" % (t1-t0))
		
	elif sys.argv[2] == 'zlib':
		sys.stderr.write("zlib module decompression test (chunksize=%d)...\n" % chunkSize)
		t0 = time.time()
		dcobj = zlib.decompressobj(zlib.MAX_WBITS | 32)
		data = sys.stdin.read(chunkSize)
		while data:
			sys.stdout.write(dcobj.decompress(data))
			data = dcobj.unconsumed_tail + sys.stdin.read(chunkSize)
		sys.stdout.write(dcobj.flush())
		t1 = time.time()
		sys.stderr.write("...completed in %f seconds\n" % (t1-t0))
		
	elif sys.argv[2] == 'gzip':
		sys.stderr.write("gzip module decompression test (chunksize=%d)...\n" % chunkSize)
		if not fileName:
			sys.stderr.write("gzip decomprssion requires filename\n")
			sys.exit(2)
		t0 = time.time()
		with gzip.GzipFile(fileName, 'rb') as f:
			data = f.read(chunkSize)
			while data:
				sys.stdout.write(data)
				data = f.read(chunkSize)
		t1 = time.time()
		sys.stderr.write("...completed in %f seconds\n" % (t1-t0))
		
elif sys.argv[1] == 'lines':
	if sys.argv[2] == 'pipe':
		sys.stderr.write("baseline (os pipe) line-based decompression test...\n")
		t0 = time.time()
		lines = 0
		proc = subprocess.Popen(["gunzip","-c"], bufsize=1, stdout=subprocess.PIPE)
		for line in iter(proc.stdout.readline, ''):
			lines += 1
			sys.stdout.write(line)
		t1 = time.time()
		sys.stderr.write("...completed %d lines in %f seconds\n" % (lines,t1-t0))
		
	elif sys.argv[2] == 'zlib':
		sys.stderr.write("zlib module decompression test (chunksize=%d)...\n" % chunkSize)
		t0 = time.time()
		lines = 0
		dcobj = zlib.decompressobj(zlib.MAX_WBITS | 32)
		data = sys.stdin.read(chunkSize)
		tail = ''
		while data:
			buf = (tail + dcobj.decompress(data)).split('\n')
			tail = buf.pop()
			for line in buf:
				lines += 1
				sys.stdout.write(line)
				sys.stdout.write('\n')
			data = dcobj.unconsumed_tail + sys.stdin.read(chunkSize)
		if tail:
			lines += 1
			sys.stdout.write(tail)
		sys.stdout.write(dcobj.flush())
		t1 = time.time()
		sys.stderr.write("...completed %d lines in in %f seconds\n" % (lines,t1-t0))
		
	elif sys.argv[2] == 'gzip': # this is the one that takes 2x as long as the other methods
		sys.stderr.write("gzip module line-based decompression test...\n")
		if not fileName:
			sys.stderr.write("gzip decomprssion requires filename\n")
			sys.exit(2)
		t0 = time.time()
		lines = 0
		with gzip.GzipFile(fileName, 'rb') as f:
			for line in f:
				lines += 1
				sys.stdout.write(line)
		t1 = time.time()
		sys.stderr.write("...completed %d lines in in %f seconds\n" % (lines,t1-t0))
