#!/usr/bin/env python

import apsw
import collections
import sys

db = apsw.Connection(sys.argv[1])
cursor = db.cursor()

# load namespace definitions
nsName = dict()
nsID = dict()
for row in cursor.execute("select namespace_id,namespace from namespace;"):
	nsName[row[0]] = row[1]
	nsID[row[1]] = row[0]

# load namespace aliases
alias = dict()
for row in cursor.execute("select namespace_id1,name1,namespace_id2,name2 from unit_name_name where namespace_id1 = namespace_id2;"):
	n1 = (row[0],row[1])
	n2 = (row[2],row[3])
	if n1 < n2:
		alias[n2] = min(n1,alias.get(n2,n1))
	elif n1 > n2:
		alias[n1] = min(n2,alias.get(n1,n2))
print "%d namespace aliases" % (len(alias),)

# load name graph
graph = collections.defaultdict(set)
for row in cursor.execute("select namespace_id1,name1,namespace_id2,name2 from unit_name_name where namespace_id1 != namespace_id2;"):
	n1 = (row[0],row[1])
	n2 = (row[2],row[3])
	while n1 in alias:
		n1 = alias[n1]
	while n2 in alias:
		n2 = alias[n2]
	if n1 != n2:
		graph[n1].add(n2)
		graph[n2].add(n1)
print "raw graph:"
print "  nodes:",len(graph)
print "  edges:",sum(len(graph[n]) for n in graph)/2


def connectedComponents(graph, only=None, skip=None, skipNode=None):
	only = set(nsID[ns] for ns in (only or nsID))
	skip = set(nsID[ns] for ns in (skip or []))
	skipNode = skipNode or set()
	nodes = set(node for node in graph if ((node[0] in only) and (node[0] not in skip)))
 	stack = list()
	while nodes:
		node = nodes.pop()
		if node in skipNode:
			continue
		component = {node}
		stack.append(node)
		while stack:
			node = stack.pop()
			for neighbor in graph[node]:
				component.add(neighbor)
				if neighbor in nodes:
					nodes.remove(neighbor)
					stack.append(neighbor)
		yield component
#connectedComponents()


def nameComponents(components):
	names = collections.defaultdict(set)
	for i,c in enumerate(components):
		for n in c:
			names[n].add(i)
	return names
#nameComponents()


def componentBreakdown(components, namespace, head=None, tail=None):
	namespaceID = nsID[namespace]
	breakdown = collections.defaultdict(int)
	for c in components:
		breakdown[sum((1 if (n[0] == namespaceID) else 0) for n in c)] += 1
	l = [("%dx%d" % (breakdown[n],n)) for n in sorted(breakdown)]
	if (head or tail) and (((head or 0) + (tail or 0)) < len(l)):
		l = l[:(head or 0)] + ["..."] + l[-(tail or 0):]
	return l
#componentBreakdown


def reverseBreakdown(names, namespace, head=None, tail=None):
	namespaceID = nsID[namespace]
	breakdown = collections.defaultdict(int)
	empty = list()
	for n in graph:
		if n[0] == namespaceID:
			breakdown[len(names.get(n,empty))] += 1
	l = [("%dx%d" % (breakdown[n],n)) for n in sorted(breakdown)]
	if (head or tail) and (((head or 0) + (tail or 0)) < len(l)):
		l = l[:(head or 0)] + ["..."] + l[-(tail or 0):]
	return l
#reverseBreakdown


def nameBreakdown(names, namespace):
	b0 = b1 = bD = bN = 0
	namespaceID = nsID[namespace]
	hits = dict()
	for n in graph:
		if n[0] == namespaceID:
			nc = names.get(n)
			if (not nc) or (len(nc) < 1):
				b0 += 1
			elif len(nc) > 1:
				bN += 1
			else:
				for c in nc:
					if c not in hits:
						b1 += 1
						hits[c] = 1
					elif hits[c] == 1:
						b1 -= 1
						bD += 2
						hits[c] = 2
					else:
						bD += 1
	return ["%dx0"%b0, "%dx1"%b1, "%dxD"%bD, "%dx2+"%bN]
#nameBreakdown


def matchBreakdown(graph, names, namespaces=None):
	namespaceIDs = set(nsID[ns] for ns in (namespaces or nsID))
	bX = b0 = b1 = bN = 0
	for n1,n2s in graph.iteritems():
		if n1[0] not in namespaceIDs:
			continue
		for n2 in n2s:
			if n2[0] not in namespaceIDs:
				continue
			if n2 <= n1:
				continue
			if (n1 not in names) or (n2 not in names):
				bX += 1
			elif len(names[n1] & names[n2]) < 1:
				b0 += 1
			elif len(names[n1] | names[n2]) > 1:
				bN += 1
			else:
				b1 += 1
	return ["%dx0"%bX, "%dx-"%b0, "%dx+"%b1, "%dx~"%bN]
#matchBreakdown()				


# generic components
components = list(connectedComponents(graph, skip={'label','description','utype'}))
names = nameComponents(components)
print "all components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# non-protein components
components = list(connectedComponents(graph, skip={'label','description','utype','ensembl_pid','refseq_pid','uniprot_pid'}))
names = nameComponents(components)
print "non-protein components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# dbID components
components = list(connectedComponents(graph, skip={'label','description','utype','ensembl_pid','refseq_pid','uniprot_pid','symbol'}))
names = nameComponents(components)
print "ID components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# ccds components
components = list(connectedComponents(graph, only={'ccds_gid'}))
names = nameComponents(components)
print "ccds components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# ensembl components
components = list(connectedComponents(graph, only={'ensembl_gid'}))
names = nameComponents(components)
print "ensembl components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# hgnc components
components = list(connectedComponents(graph, only={'hgnc_gid'}))
names = nameComponents(components)
print "hgnc components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# ncbi/entrez components
components = list(connectedComponents(graph, only={'entrez_gid'}))
names = nameComponents(components)
print "ncbi/entrez components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# N,E components
skipEnsembl = set()
components1 = list(connectedComponents(graph, only={'entrez_gid'}))
for c in components1:
	for n in c:
		if n[0] == nsID['ensembl_gid']:
			skipEnsembl.add(n)
components2 = list(connectedComponents(graph, only={'ensembl_gid'}, skipNode=skipEnsembl))
components = components1 + components2
components1 = components2 = None
names = nameComponents(components)
print "N,E sequential components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# C,H,N,E components
skipHGNC = set()
skipEntrez = set()
skipEnsembl = set()
components1 = list(connectedComponents(graph, only={'ccds_gid'}))
for c in components1:
	for n in c:
		if n[0] == nsID['ensembl_gid']:
			skipEnsembl.add(n)
		elif n[0] == nsID['hgnc_gid']:
			skipHGNC.add(n)
		elif n[0] == nsID['entrez_gid']:
			skipEntrez.add(n)
components2 = list(connectedComponents(graph, only={'hgnc_gid'}, skipNode=skipHGNC))
for c in components2:
	for n in c:
		if n[0] == nsID['ensembl_gid']:
			skipEnsembl.add(n)
		elif n[0] == nsID['entrez_gid']:
			skipEntrez.add(n)
components3 = list(connectedComponents(graph, only={'entrez_gid'}, skipNode=skipEntrez))
for c in components3:
	for n in c:
		if n[0] == nsID['ensembl_gid']:
			skipEnsembl.add(n)
components4 = list(connectedComponents(graph, only={'ensembl_gid'}, skipNode=skipEnsembl))
components = components1 + components2 + components3 + components4
components1 = components2 = components3 = components4 = None
names = nameComponents(components)
print "C,H,N,E sequential components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# E,N,H,C components
skipEntrez = set()
skipHGNC = set()
skipCCDS = set()
components1 = list(connectedComponents(graph, only={'ensembl_gid'}))
for c in components1:
	for n in c:
		if n[0] == nsID['ccds_gid']:
			skipCCDS.add(n)
		elif n[0] == nsID['hgnc_gid']:
			skipHGNC.add(n)
		elif n[0] == nsID['entrez_gid']:
			skipEntrez.add(n)
components2 = list(connectedComponents(graph, only={'entrez_gid'}, skipNode=skipEntrez))
for c in components2:
	for n in c:
		if n[0] == nsID['ccds_gid']:
			skipCCDS.add(n)
		elif n[0] == nsID['hgnc_gid']:
			skipHGNC.add(n)
components3 = list(connectedComponents(graph, only={'hgnc_gid'}, skipNode=skipHGNC))
for c in components3:
	for n in c:
		if n[0] == nsID['ccds_gid']:
			skipCCDS.add(n)
components4 = list(connectedComponents(graph, only={'ccds_gid'}, skipNode=skipCCDS))
components = components1 + components2 + components3 + components4
components1 = components2 = components3 = components4 = None
names = nameComponents(components)
print "E,N,H,C sequential components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# analyze E+N simultaneous component breakdown
components = list(connectedComponents(graph, only={'ensembl_gid','entrez_gid'}))
names = nameComponents(components)
print "E+N simultaneous components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

# analyze C+E+H+N simultaneous component breakdown
components = list(connectedComponents(graph, only={'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
names = nameComponents(components)
print "C+E+H+N simultaneous components:",len(components)
print "   edge:",", ".join(matchBreakdown(graph,names,{'ccds_gid','ensembl_gid','entrez_gid','hgnc_gid'}))
print " > labl:",", ".join(componentBreakdown(components,'label',4,4))
print " < ccds:",", ".join(nameBreakdown(names,'ccds_gid'))
print " < ensg:",", ".join(nameBreakdown(names,'ensembl_gid'))
print " < hgnc:",", ".join(nameBreakdown(names,'hgnc_gid'))
print " < ncbi:",", ".join(nameBreakdown(names,'entrez_gid'))

