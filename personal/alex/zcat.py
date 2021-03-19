#!/usr/bin/env python

import zlib

def zcat(file_name, split="\n", chunk=1*1024*1024):
	dc = zlib.decompressobj(32+zlib.MAX_WBITS) # autodetect gzip or zlib header
	with open(file_name,'rb') as f:
		text = ""
		loop = True
		while loop:
			data = f.read(chunk)
			if data:
				text += dc.decompress(data)
				data = None
			else:
				text += dc.flush()
				loop = False
			if text:
				lines = text.split(split)
				i,x = 0,len(lines)-1
				text = lines[x]
				while i < x:
					yield lines[i]
					i += 1
				lines = None
		#while data remains
		if text:
			yield text
#zcat()


if __name__ == "__main__":
	import sys
	for line in zcat(sys.argv[1],"\n",int(sys.argv[2])):
		print line
