#!/usr/bin/env python

import time
import sys

if sys.argv[1] == 'c':
	t0 = time.time()
	s = ""
	d = " "
	for i in xrange(0,10000000):
		s += ("%d%s" % (i,d))
	s = s[:-len(d)]
	t1 = time.time()
	print "string concat: %1.1fs" % (t1-t0)
elif sys.argv[1] == 'j':
	t0 = time.time()
	l = list()
	d = " "
	for i in xrange(0,10000000):
		l.append("%d" % i)
	s = d.join(l)
	l = None
	t1 = time.time()
	print "list join: %1.1fs" % (t1-t0)
else:
	print "c or j"
