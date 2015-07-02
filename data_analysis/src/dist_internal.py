#!/usr/bin/python

import sys

if len(sys.argv) != 2:
	print "must provide exactly 1 .data file"
	exit(1)

def dist(m1, m2):
	count = 0
	for i in range(len(m1)):
		if m1[i] != m2[i] and m1[i] != 'N' and m2[i] != 'N':
			count += 1
	return count

mlen = 0
d = []
with open(sys.argv[1], "r") as f:
	for line in f:
		if line[0] == "#" or line[0] == " ":
			continue
		
		motif = line.split()[0]
		if mlen == 0:
			mlen = len(motif)
		d.append(motif)

if len(d) == 0:
	print "-"
	exit(0)
elif len(d) == 1:
	print "0"
	exit(0)

tot = 0
n   = 0
for i in range(len(d)):
	for j in range(i+1, len(d)):
		tot += dist(d[i], d[j])
		n += 1

print float(tot) / (n*mlen)
