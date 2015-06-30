#!/usr/bin/python

import sys

if len(sys.argv) < 2:
	print "must provide at least 1 .data file"
	exit(1)

for filename in sys.argv[1:]:
	GC  = 0
	CpG = 0
	tot = 0
	with open(filename, "r") as f:
		for line in f:
			if line[0] == '#' or line[0] == ' ':
				continue

			motif = line.split()[0]

			prev = 'X'
			for c in motif:
				if c == 'C' or c == 'G':
					GC += 1
				if prev == 'C' and c == 'G':
					CpG += 1
				
				tot += 1
				prev = c

	print filename
	print "tot:", tot, "GC%", (float(GC) / tot) * 100.0, "CpG%", (float(CpG)*2.0 / tot) * 100.0