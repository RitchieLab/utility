#!/usr/bin/env python

import argparse
import collections
import csv
import gzip
import os
import sys


def topoSort(neighbors, nodes=None):
	visited = set()
	empty = tuple()
	for node in (neighbors or nodes):
		for nextNode in _topoSort_recurse(node, neighbors, visited, empty):
			yield nextNode
#topoSort()


def _topoSort_recurse(node, neighbors, visited, empty):
	if node not in visited:
		visited.add(node)
		for neighbor in neighbors.get(node, empty):
			for nextNode in _topoSort_recurse(neighbor, neighbors, visited, empty):
				yield nextNode
		yield node
#_topoSort_recurse()


def colcat(
		inputPaths,
		uidPaths,
		uidColumns,
		repeatUIDs,
		missingCode,
		outputPath,
		inDelimChar,
		outDelimChar,
		inQuoteChar,
		outQuoteChar,
		inEscapeChar,
		outEscapeChar,
		inTermChar,
		outTermChar,
		compressed,
		outCompressed,
		force,
		verbose
):
	msg = (sys.stderr if (outputPath or '-') == '-' else sys.stdout)
	columnUID = { uidColumns[u]:u for u in xrange(len(uidColumns or list())) }
	
	if verbose:
		msg.write("input files: '%s'\n" % ("', '".join(inputPaths)))
		if not uidColumns:
			msg.write("UID columns: N/A\n")
		else:
			if uidPaths:
				msg.write("UID files: '%s'\n" % ("', '".join(uidPaths)))
			else:
				msg.write("UIDs from input files\n")
			msg.write("UID columns: %s\n" % (", ".join(("%s" % c) for c in uidColumns)))
			msg.write("repeat UID columns: %s\n" % ("yes" if repeatUIDs else "no"))
		msg.write("fill in missing data: %s\n" % (("yes (%s)" % missingCode) if missingCode != None else "no"))
		msg.write("compressed input: %s\n" % (("yes" if compressed == '' else "optional (%s)" % compressed) if compressed != None else "no"))
		msg.write("compressed output: %s\n" % ("yes" if outCompressed else "no"))
	
	# open CSV readers for all input files
	if verbose:
		msg.write("opening input%s files ...\n" % (" and UID" if (uidColumns and uidPaths) else ""))
	inputs = range(0,len(inputPaths))
	inputFiles = list()
	inputCSVs = list()
	inputQueue = list()
	inputCols = list()
	inputRows = list()
	inputSkip = list()
	inputFill = list()
	uidFiles = list()
	uidCSVs = list()
	uidQueue = list()
	for i in inputs:
		if (inputPaths[i] or '-') == '-':
			inputFiles.append(gzip.open(sys.stdin,'rb') if compressed else sys.stdin)
		elif (compressed != None) and (os.path.exists(inputPaths[i]+compressed)):
			inputFiles.append(gzip.open(inputPaths[i]+compressed,'rb'))
		elif os.path.exists(inputPaths[f]):
			inputFiles.append(open(inputPaths[f],'rb'))
		else:
			exit("ERROR: could not find input file #%d (%s)" % (i,inputPaths[i]))
		inputCSVs.append(csv.reader(
				inputFiles[i],
				delimiter=inDelimChar,
				escapechar=inEscapeChar,
				quotechar=inQuoteChar,
				lineterminator=inTermChar
		))
		inputQueue.append(collections.deque())
		inputCols.append(0)
		inputRows.append(0)
		inputSkip.append(list())
		inputFill.append(list())
		
		if uidColumns and uidPaths:
			if (uidPaths[i] or '-') == '-':
				uidFiles.append(gzip.open(sys.stdin,'rb') if compressed else sys.stdin)
			elif (compressed != None) and (os.path.exists(uidPaths[i]+compressed)):
				uidFiles.append(gzip.open(uidPaths[i]+compressed,'rb'))
			elif os.path.exists(uidPaths[f]):
				uidFiles.append(open(uidPaths[f],'rb'))
			else:
				exit("ERROR: could not find UID file #%d (%s)" % (i,uidPaths[i]))
			uidCSVs.append(csv.reader(
				uidFiles[i],
				delimiter=inDelimChar,
				escapechar=inEscapeChar,
				quotechar=inQuoteChar,
				lineterminator=inTermChar
			))
		uidQueue.append(collections.deque())
	#foreach input
	if verbose:
		msg.write("... OK\n")
	
	# open CSV writer for output file
	if verbose:
		msg.write("concatenating columns...\n")
	if (outputPath or '-') == '-':
		outputFile = gzip.GzipFile(None, 'wb', compresslevel=6, fileobj=sys.stdout) if outCompressed else sys.stdout
	elif (not force) and os.path.exists(outputPath):
		exit("ERROR: output file '%s' already exists" % outputPath)
	elif outCompressed:
		outputFile = gzip.open(outputPath, 'wb', compresslevel=6)
	else:
		outputFile = open(outputPath, 'wb')
	with outputFile:
		outputCSV = csv.writer(
				outputFile,
				delimiter=outDelimChar,
				escapechar=outEscapeChar,
				quotechar=outQuoteChar,
				lineterminator=outTermChar
		)
		
		# search for matching UIDs until no lines remain
		uidFound = collections.defaultdict(set)
		mergeQueue = collections.deque()
		remain = True
		while remain:
			minQueue = min(len(queue) for queue in inputQueue)
			assert(minQueue == 0)
			if uidColumns:
				# read lines from each input until a UID match is found
				while remain and not mergeQueue:
					remain = False
					for i in inputs:
						if len(inputQueue[i]) == minQueue:
							try:
								data = inputCSVs[i].next()
								if uidPaths:
									uidData = uidCSVs[i].next()
									uid = tuple(uidData[c] for c in uidColumns)
								else:
									uid = tuple(data[c] for c in uidColumns)
							except StopIteration:
								continue
							remain = True
							if not inputCols[i]:
								inputCols[i] = len(data)
							inputRows[i] += 1
							inputQueue[i].append(data)
							uidQueue[i].append(uid)
							uidFound[uid].add(i)
							if len(uidFound[uid]) == len(inputs):
								mergeQueue.append(uid)
								break
						#if input queue isn't ahead of the others
					#foreach input
					if remain:
						minQueue += 1
				#while lines remain and no match
			else:
				assert(max(len(queue) for queue in inputQueue) == 0)
				# read one line from each input
				remain = False
				eof = False
				for i in inputs:
					try:
						data = inputCSVs[i].next()
					except StopIteration:
						eof = True
						continue
					remain = True
					if not inputCols[i]:
						inputCols[i] = len(data)
					inputRows[i] += 1
					inputQueue[i].append(data)
					uidQueue[i].append(0)
					uidFound[0].add(i)
				#foreach input
				if remain:
					minQueue += 1
				if not eof:
					mergeQueue.append(0)
			#if uidColumns
			assert(len(mergeQueue) <= 1)
			
			# if there are any lines before the full match or EOF, decide what to do with them
			if minQueue > len(mergeQueue):
				if missingCode:
					if uidColumns:
						# gather UID ordering contraints from each input
						uidAfter = dict()
						uidList = list()
						for i in inputs:
							prev = None
							for uid in uidQueue[i]:
								if mergeQueue and uid == mergeQueue[-1]:
									break
								if uid not in uidAfter:
									uidAfter[uid] = set()
									uidList.append(uid)
								if prev:
									uidAfter[uid].add(prev)
								prev = uid
							#foreach uid before match
						#foreach input
						
						# queue partial matches in topological-sorted order
						try:
							uid = mergeQueue.popleft()
						except IndexError:
							uid = None
						mergeQueue.extend(topoSort(uidAfter, uidList))
						if uid:
							mergeQueue.append(uid)
					else:
						mergeQueue.append(0)
					#if uidColumns
				else:
					# with no missing code, just discard the extras
					for i in inputs:
						while (len(uidQueue[i]) > 0) and ((not uidColumns) or (not mergeQueue) or (uidQueue[i][0] != mergeQueue[-1])):
							inputQueue[i].popleft()
							uid = uidQueue[i].popleft()
							uidFound[uid].remove(i)
							inputSkip[i].append(uid)
							msg.write("WARNING: dropped unmatched line #%d from input #%d (%s)%s\n" % (
								inputRows[i],i+1,inputPaths[i],((" with ID %s" % uid) if uidColumns else "")
							))
					#foreach input
				#if missingCode
			#if minQueue > len(mergeQueue)
			
			# merge rows in order
			while mergeQueue:
				uid = mergeQueue.popleft()
				merged = list()
				for i in inputs:
					if i in uidFound[uid]:
						# use the real data if available from this input
						data = inputQueue[i].popleft()
						assert(uid == uidQueue[i][0])
						uidQueue[i].popleft()
						uidFound[uid].remove(i)
						if not uidColumns:
							merged.extend(data)
						elif i == 0 or repeatUIDs:
							if uidPaths:
								merged.extend(uid)
							merged.extend(data)
						else:
							merged.extend(data[c] for c in xrange(len(data)) if c not in columnUID)
					else:
						# fill in data if missing from the input
						assert((not uidQueue[i]) or (uid != uidQueue[i][0]))
						if not uidColumns:
							merged.extend(missingCode for c in xrange(inputCols[i]))
						elif i == 0 or repeatUIDs:
							if uidPaths:
								merged.extend(uid)
								merged.extend(missingCode for c in xrange(inputCols[i]))
							else:
								merged.extend((uid[columnUID[c]] if c in columnUID else missingCode) for c in xrange(inputCols[i]))
						else:
							merged.extend(missingCode for c in xrange(inputCols[i]) if c not in columnUID)
						inputFill[i].append(uid)
						msg.write("WARNING: filled in data missing from input #%d (%s)%s\n" % (
							i+1,inputPaths[i],((" for ID %s" % uid) if uidColumns else "")
						))
					#if input in uidFound[uid]
				#foreach input
				assert(not uidFound[uid])
				del uidFound[uid]
				outputCSV.writerow(merged)
			#while mergeQueue
		#while lines remain
	#with outputFile
	if verbose:
		msg.write("... OK\n")
	
	# close all input files
	for file in inputFiles:
		file.close()
#colcat()


if __name__ == "__main__":
	versName,versMaj,versMin,versRev,versDate = sys.argv[0],0,1,0,'2012-04-19'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		#formatter_class=argparse.RawDescriptionHelpFormatter,
		description="%s version %s" % (versName, version),
		epilog="""
example: %(prog)s TODO
"""
	)
	parser.add_argument('inputs', type=str, metavar='file', action='store', nargs='+', default=['-'],
		help="input filename(s); may contain wildcards to be expanded using the specified sequence"
	)
	parser.add_argument('-u', '--uid-files', type=str, metavar='file', action='store', nargs='+',
		help="filename(s) which contain unique IDs for the respective input file(s); "+
		"may contain wildcards to be expanded using the specified sequence"
	)
	parser.add_argument('-w', '--wildcard', type=str, metavar='char', action='store', default='#',
		help="wildcard character(s) to expand in input or UID filename(s) (default: #)"
	)
	parser.add_argument('-s', '--sequence', type=str, metavar='list/range', action='store',
		help="sequence with which to expand input or UID filename wildcards (default: no wildcard expansion)"
	)
	parser.add_argument('-c', '--uid-columns', type=str, metavar='list/range', action='store',
		help="column(s) containing unique ID values with which to identify corresponding lines of each input; "+
		"taken from --uid-files if specified, otherwise from the inputs"
	)
	parser.add_argument('-r', '--repeat-uids', action='store_true',
		help="include the matching unique ID column(s) from every input "+
		"(default: skip ID column(s) from inputs after the first)"
	)
	parser.add_argument('-m', '--missing-code', type=str, metavar='char(s)', action='store',
		help="missing-data code with which to fill in for lines missing in some inputs "+
		"(default: skip lines without matches in all inputs)"
	)
	parser.add_argument('-o', '--output', type=str, metavar='file', action='store', default='-',
		help="output filename (default: stdout)"
	)
	parser.add_argument('-d', '--delimiter', type=str, metavar='char', action='store', default='\t',
		help="character which delimits input fields (default: tab)"
	)
	parser.add_argument('-D', '--output-delimiter', type=str, metavar='char', action='store', default=None,
		help="character with which to delimit output fields (default: same as input delimiter)"
	)
	parser.add_argument('-q', '--quote', type=str, metavar='char', action='store', default=None,
		help="character which optionally encloses input fields (default: none)"
	)
	parser.add_argument('-Q', '--output-quote', type=str, metavar='char', action='store', default=None,
		help="character with which to optionally enclose output fields (default: same as input quote)"
	)
	parser.add_argument('-e', '--escape', type=str, metavar='char', action='store', default=None,
		help="character which escapes delimiters and quotes in input data (default: none)"
	)
	parser.add_argument('-E', '--output-escape', type=str, metavar='char', action='store', default=None,
		help="character with which to escape delimiters and quotes in output data (default: same as input escape)"
	)
	parser.add_argument('-t', '--terminator', type=str, metavar='char', action='store', default='\n',
		help="character(s) which terminate input lines (default: \\n)"
	)
	parser.add_argument('-T', '--output-terminator', type=str, metavar='char', action='store', default=None,
		help="character(s) with which to terminate output lines (default: same as input terminator)"
	)
	parser.add_argument('-z', '--compressed', type=str, metavar='ext', nargs='?', action='store', default=None, const='',
		help="read compressed input data; with an optional extension, look for input filename(s) plus extension "+
		"from which to read compressed data if available, otherwise use base filename(s) as normal"
	)
	parser.add_argument('-Z', '--output-compressed', action='store_true',
		help="write compressed output data"
	)
	parser.add_argument('-f', '--force', action='store_true',
		help="force overwrite of existing output files"
	)
	parser.add_argument('-v', '--verbose', action='store_true',
		help="print verbose progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# parse and apply wildcard sequence, if specified
	if args.sequence:
		sequence = list()
		for substitute in args.sequence.split(','):
			# try interpreting it as a range
			dash = substitute.find('-',1)
			colon = substitute.find(':',dash+1)
			colon = colon if colon > dash else None
			try:
				start = int(substitute[0:dash])
				stop = int(substitute[dash+1:colon])+1
				step = int(substitute[colon+1:]) if colon else 1
				sequence.extend(str(n) for n in xrange(start,stop,step))
			except ValueError:
				# if int() complains then it's not a numeric range; use it as-is
				sequence.append(substitute)
		#foreach substitute
		
		# do wildcard expansions backwards so the expanding list doesn't break the iteration
		for n in xrange(len(args.inputs)-1,-1,-1):
			if args.wildcard in args.inputs[n]:
				expansion = list()
				for substitute in sequence:
					expansion.append(args.inputs[n].replace(args.wildcard, substitute))
				args.inputs[n:n+1] = expansion
		
		if args.uid_files:
			for n in xrange(len(args.uid_files)-1,-1,-1):
				if args.wildcard in args.uid_files[n]:
					expansion = list()
					for substitute in sequence:
						expansion.append(args.uid_files[n].replace(args.wildcard, substitute))
					args.uid_files[n:n+1] = expansion
	#if args.sequence
	
	# parse ID columns, if specified
	uidColumns = None
	if args.uid_columns:
		uidColumns = list()
		for column in args.uid_columns.split(','):
			# try interpreting it as a range
			dash = column.find('-',1)
			colon = column.find(':',dash+1)
			colon = colon if colon > dash else None
			try:
				start = int(column[0:dash])
				stop = int(column[dash+1:colon])+1
				step = int(column[colon+1:]) if colon else 1
				uidColumns.extend(range(start-1,stop-1,step))
			except ValueError:
				# if int() complains on the range components then it's not a numeric range; use it as-is
				try:
					uidColumns.append(int(column)-1)
					if uidColumns[-1] < 0:
						raise ValueError()
				except ValueError:
					sys.exit("ERROR: invalid ID column index '%s', must be a positive integer" % column)
		#foreach column
	#if args.columns
	
	# validate ID file count, if specified
	if args.uid_files:
		if len(args.uid_files) != len(args.inputs):
			sys.exit("ERROR: the number of UID files (%d) does not match the number of input files (%d)" % (len(args.uid_files),len(args.inputs)))
	
	# run
	colcat(
		args.inputs,
		args.uid_files,
		uidColumns,
		args.repeat_uids,
		args.missing_code,
		args.output,
		args.delimiter,
		args.output_delimiter or args.delimiter,
		args.quote,
		args.output_quote or args.quote,
		args.escape,
		args.output_escape or args.escape,
		args.terminator,
		args.output_terminator or args.terminator or '\n',
		args.compressed,
		args.output_compressed,
		args.force,
		args.verbose
	)
#__main__
