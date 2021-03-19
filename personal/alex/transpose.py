#!/usr/bin/env python

import argparse
import csv
import gzip
import math
import os
import sys


def transpose(
		inputPath, outputPath,
		numHeaders, numLabels, numFields,
		mergeHeaders, repeatLabels,
		delimChar, quoteChar, escChar, termChar,
		delimCharOut, quoteCharOut, escCharOut, termCharOut,
		extraDelim, memLimit, compressed, compressOut,
		verbose
):
	msg = (sys.stderr if (outputPath or '-') == '-' else sys.stdout)
	
	# calculate object overhead on this system
	memList = sys.getsizeof([])
	memStr = sys.getsizeof('')
	
	if verbose:
		import time
		memLimitStr = "no"
		for unit in (('TB',1024*1024*1024*1024), ('GB',1024*1024*1024), ('MB',1024*1024), ('KB',1024), ('B',1)):
			if memLimit >= unit[1] * 0.9:
				memLimitStr = "%1.1f %s" % (1.0*memLimit/unit[1],unit[0])
				break
		msg.write("processing input file (%s buffer limit):\n" % memLimitStr)
		t0 = time.time()
	with (sys.stdin if (inputPath or '-') == '-' else open(inputPath,'rb')) as inputFile:
		with (sys.stdout if (outputPath or '-') == '-' else open(outputPath,'wb')) as outputFile:
			out = csv.writer(outputFile, delimiter=delimCharOut, escapechar=escCharOut, quotechar=quoteCharOut, lineterminator=termCharOut)
			# initialize data
			inputSize = None
			fieldRange = range(0,numFields)
			iColNum = None
			iCol0 = 0
			iCol1 = None
			oRowNum = None
			oRow0 = 0
			oRow1 = None
			# keep doing passes until the last column is transposed
			p = 0
			while True:
				p += 1
				if verbose:
					msg.write("  pass %d%s ..." % (p,(" of ~%d" % (p + int(math.ceil(1.0 * (iColNum - iCol0) / (iCol1 - iCol0))) - 1)) if p > 1 else ""))
					msg.flush()
					t0p = time.time()
				# reset the input file, if possible
				try:
					if verbose and (p == 1):
						inputFile.seek(0, os.SEEK_END)
						inputSize = inputFile.tell()
					inputFile.seek(0, os.SEEK_SET)
					canSeek = True
				except IOError as e:
					if e.errno != 29:
						raise e
					canSeek = False
				decile = 0
				memUsage = 0
				oRows = None
				# read and transpose selected columns from each line
				iRow = 0
				for iCols in csv.reader(inputFile, delimiter=delimChar, escapechar=escChar, quotechar=quoteChar, lineterminator=termChar):
					if extraDelim and iCols[-1] == "":
						iCols.pop()
					
					# during the first pass, audit the input data
					if p == 1:
						if iRow == 0:
							iCol1 = iColNum = len(iCols)
							if (iColNum - numLabels) % numFields != 0:
								msg.write(" ERROR: found %d input columns, expected %d+%dx\n" % (iColNum,numLabels,numFields))
								sys.exit(1)
							oRow1 = oRowNum = numLabels + (iColNum - numLabels) / numFields
						elif len(iCols) < iColNum:
							msg.write(" ERROR: found %d input columns on line %d, expected at least %d\n" % (len(iCols),iRow+1,iCol1))
							sys.exit(1)
					
					# on the first line, reset the output buffer
					if iRow == 0:
						oRows = [ [] for r in xrange(0,oRow1-oRow0) ]
						memUsage = memList * (1 + len(oRows))
					
					# transpose columns
					iCol = iCol0
					oRow = oRow0
					while iCol < iCol1:
						if memLimit:
							if iCol < numLabels:
								memUsage += memStr*numFields + len(iCols[iCol])
							else:
								memUsage += memStr*numFields + sum(len(iCols[iCol+f]) for f in fieldRange)
							
							while memUsage > memLimit:
								if not canSeek:
									msg.write(" ERROR: memory limit exceeded, input stream cannot be reset for another pass\n")
									sys.exit(1)
								iCol1 -= (1 if iCol1 <= numLabels else numFields)
								if iCol1 <= iCol0:
									msg.write(" ERROR: memory limit exceeded on single column #%d\n" % (iCol0+1))
									sys.exit(1)
								oRow1 -= 1
								memUsage -= memList + memStr*len(oRows[oRow1-oRow0]) + sum(len(s) for s in oRows[oRow1-oRow0])
								del oRows[oRow1-oRow0]
							
							if iCol >= iCol1:
								break
						#if memLimit
						
						if iCol < numLabels:
							oRows[oRow-oRow0].extend(("" if f else iCols[iCol]) for f in fieldRange)
							iCol += 1
						elif iRow < numHeaders:
							#TODO
							oRows[oRow-oRow0].extend(iCols[iCol+f] for f in fieldRange)
							iCol += numFields
						else:
							oRows[oRow-oRow0].extend(iCols[iCol+f] for f in fieldRange)
							iCol += numFields
						oRow += 1
					#while iCol < iCol1
					iRow += 1
					
					if inputSize:
						progress = int(10.0 * inputFile.tell() / inputSize) * 10
						if progress > decile:
							decile = progress
							if progress >= 100:
								msg.write(" ...")
							else:
								msg.write(" %d%%" % decile)
							msg.flush()
				#while inputFile has lines
				for r in xrange(0,oRow1-oRow0):
					out.writerow(oRows[r])
				if verbose:
					t1p = time.time()
					msg.write(" OK: %d columns (#%d-%d) in %1.1fs\n" % (iCol1-iCol0,iCol0+1,iCol1,t1p-t0p))
					msg.flush()
				if iCol1 >= iColNum:
					break
				iCol0,iCol1 = iCol1,iCol1+iCol1-iCol0
				if iCol1 >= numLabels:
					iCol1 *= 1.0 * memLimit / memUsage
					iCol1 = int(math.ceil(1.0 * (iCol1 - numLabels) / numFields)) * numFields + numLabels
				iCol1 = min(iCol1,iColNum)
				oRow0,oRow1 = oRow1,(iCol1 if iCol1 < numLabels else (numLabels + (iCol1 - numLabels) / numFields))
			#while columns remain to transpose
		#with out
	#with inputFile
	if verbose:
		t1 = time.time()
		msg.write("transposition complete in %1.1fs\n" % (t1-t0))
#transpose()


if __name__ == "__main__":
	versName,versMaj,versMin,versRev,versDate = 'transpose',0,1,0,'2012-04-04'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		#formatter_class=argparse.RawDescriptionHelpFormatter,
		description="%s version %s" % (versName, version),
		epilog="""
example: %(prog)s input.txt -o output.txt -r 1 -c 2 -f 3 -d " " -v
"""
	)
	parser.add_argument('input', type=str, metavar='file', action='store', default='-',
		help="input filename (default: stdin)"
	)
	parser.add_argument('-o', '--output', type=str, metavar='file', action='store', default='-',
		help="output filename (default: stdout)"
	)
	parser.add_argument('-r', '--header-rows', type=int, metavar='num', action='store', default=0,
		help="number of header rows (default: 0)"
	)
	parser.add_argument('-c', '--label-cols', type=int, metavar='num', action='store', default=0,
		help="number of label columns (default: 0)"
	)
	parser.add_argument('-f', '--fields', type=int, metavar='num', action='store', default=1,
		help="number of fields to be transposed as one cell (default: 1)"
	)
	parser.add_argument('-g', '--merge-headers', action='store_true',
		help="when --header-rows > 0 and --fields > 1, merge empty or repeated header fields (default: all header fields transposed into label columns)"
	)
	parser.add_argument('-p', '--repeat-labels', action='store_true',
		help="when --label-cols > 0 and --fields > 1, generate filler fields in header rows by repeating the label(s) (default: filler fields are blank)"
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
	parser.add_argument('-x', '--extra-delimiter', action='store_true',
		help="ignore one trailing delimiter on each input line (default: no)"
	)
	parser.add_argument('-m', '--memory-limit', type=str, metavar='num', action='store', default=None,
		help="target size for transpose buffer, with optional unit suffix: b,k,m,g,t (default: no limit)"
	)
	parser.add_argument('-z', '--compressed', type=str, metavar='ext', nargs='?', action='store', default=None, const='',
		help="read compressed input data; with an optional extension, look for input filename plus extension "+
		"from which to read compressed data if available, otherwise use base filename as normal"
	)
	parser.add_argument('-Z', '--output-compressed', action='store_true',
		help="write compressed output data"
	)
	parser.add_argument('-v', '--verbose', action='store_true',
		help="print verbose progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# run
	if args.memory_limit:
		m = args.memory_limit.strip()
		if m[-1] == 'b' or m[-1] == 'B':
			m = m[:-1]
		if m[-1] == 'k':
			m = float(m[:-1]) * (1000)
		elif m[-1] == 'K':
			m = float(m[:-1]) * (1024)
		elif m[-1] == 'm':
			m = float(m[:-1]) * (1000 * 1000)
		elif m[-1] == 'M':
			m = float(m[:-1]) * (1024 * 1024)
		elif m[-1] == 'g':
			m = float(m[:-1]) * (1000 * 1000 * 1000)
		elif m[-1] == 'G':
			m = float(m[:-1]) * (1024 * 1024 * 1024)
		elif m[-1] == 't':
			m = float(m[:-1]) * (1000 * 1000 * 1000 * 1000)
		elif m[-1] == 'T':
			m = float(m[:-1]) * (1024 * 1024 * 1024 * 1024)
		args.memory_limit = long(m)
	transpose(
			args.input,
			args.output,
			args.header_rows,
			args.label_cols,
			args.fields,
			args.merge_headers,
			args.repeat_labels,
			args.delimiter,
			args.quote,
			args.escape,
			args.terminator,
			args.output_delimiter or args.delimiter,
			args.output_quote or args.quote,
			args.output_escape or args.escape,
			args.output_terminator or args.terminator or '\n',
			args.extra_delimiter,
			args.memory_limit,
			args.compressed,
			args.output_compressed,
			args.verbose
	)
#__main__
