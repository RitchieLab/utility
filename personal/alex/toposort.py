#!/usr/bin/env python

def topoSort(neighbors, nodes=None):
	visited = set()
	empty = tuple()
	for node in (neighbors or nodes):
		for nextNode in _topoSort_recurse(node, neighbors, visited, empty):
			yield nextNode


def _topoSort_recurse(node, neighbors, visited, empty):
	if node not in visited:
		visited.add(node)
		for neighbor in neighbors.get(node, empty):
			for nextNode in _topoSort_recurse(neighbor, neighbors, visited, empty):
				yield nextNode
		yield node


import sys
import time

before = dict()
print "sequences:"
for seq in sys.argv[1:]:
	print "  ",seq
	for i in xrange(0,len(seq)):
		if seq[i] not in before:
			before[seq[i]] = set()
		if 0:
			before[seq[i]] |= set(seq[0:i])
		elif i > 0:
			before[seq[i]].add(seq[i-1:i])

print "requirements:"
for item in before:
	print "  ",item,": ",before[item]

print "ordering:"
print "  ",tuple(topoSort(before))

n = 10000000
print "timing",n,"runs:"
t0 = time.time()
while n > 0:
	topoSort(before)
	n -= 1
t1 = time.time()
print "  ",t1-t0,"s"
