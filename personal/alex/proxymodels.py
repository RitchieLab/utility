#!/usr/bin/env python

import collections
import sys

args = sys.argv
only = label = False
for a in xrange(3,len(args)):
	if args[a].lower() == 'only':
		only = True
	elif args[a].lower() == 'label':
		label = True
	else:
		args = list()
		break
if len(args) < 3:
	print "usage: %s proxyfile modelfile [only] [label] > outputfile" % args[0]
	sys.exit(2)

proxy = collections.defaultdict(set)
with (sys.stdin if args[1] == "-" else open(args[1],'rU')) as proxyFile:
	for line in proxyFile:
		words = list(w.strip() for w in line.split(None))
		if len(words) >= 2:
			proxy[words[0]].add(words[1])

blank = ['']
n = 0
with (sys.stdin if args[2] == "-" else open(args[2],'rU')) as modelFile:
	for line in modelFile:
		words = list(w.strip() for w in line.split(None))
		if len(words) >= 2:
			n += 1
			if not only:
				proxy[words[0]].add(words[0])
				proxy[words[1]].add(words[1])
			for p0 in (proxy[words[0]] or blank):
				for p1 in (proxy[words[1]] or blank):
					if label:
						print p0,p1,("orig%d" if (p0==words[0]) and (p1==words[1]) else "prox%d") % n
					else:
						print p0,p1
