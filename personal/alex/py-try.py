#!/usr/bin/env python

# results: try/except is faster than explicit in-dict or in-list testing
#          as long as the exception rate is no more than about 10%


import sys
import time

rate = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
size = int(sys.argv[2]) if len(sys.argv) > 2 else 1000000
testDict = { n:1 for n in xrange(size) if (n % (1/rate)) < ((n+1) % (1/rate)) }
testList = [ 1 for n in testDict ]

print "testing %1.1f%% exception rate over %d iterations" % (100*(1.0-(1.0*len(testList)/size)),size)

ok = 0
bad = 0
t0 = time.time()
for n in xrange(size):
	if n in testDict:
		ok += testDict[n]
	else:
		bad += 1
t1 = time.time()
print "if in dict: %d of %d hits, %d of %d misses in %1.1f seconds" % (ok,len(testDict),bad,size-len(testDict),t1-t0)

ok = 0
bad = 0
t0 = time.time()
for n in xrange(size):
	try:
		ok += testDict[n]
	except KeyError:
		bad += 1
t1 = time.time()
print "try dict[key]: %d of %d hits, %d of %d misses in %1.1f seconds" % (ok,len(testDict),bad,size-len(testDict),t1-t0)

ok = 0
bad = 0
t0 = time.time()
for n in xrange(size):
	if n < len(testList):
		ok += testList[n]
	else:
		bad += 1
t1 = time.time()
print "if n < len(list): %d of %d hits, %d of %d misses in %1.1f seconds" % (ok,len(testList),bad,size-len(testList),t1-t0)

ok = 0
bad = 0
t0 = time.time()
for n in xrange(size):
	try:
		ok += testList[n]
	except IndexError:
		bad += 1
t1 = time.time()
print "try list[n]: %d of %d hits, %d of %d misses in %1.1f seconds" % (ok,len(testList),bad,size-len(testList),t1-t0)
