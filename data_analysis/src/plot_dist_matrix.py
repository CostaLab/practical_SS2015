#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import pylab

if len(sys.argv) < 2:
	print "must provide at least 1 .data file"
	exit(1)

def dist(m1, m2):
	count = 0
	for i in range(len(m1)):
		if m1[i] != m2[i] and m1[i] != 'N' and m2[i] != 'N':
			count += 1
	return count

for filename in sys.argv[1:]:
	print "# plotting distance matrix for", filename
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
		continue

	m = np.empty((len(d), len(d)), dtype=float)
	for i in range(len(d)):
		for j in range(len(d)):
			if i == j:
				m[i][j] = 0.0
			elif j < i:
				# we've already calculated it, let's copy
				m[i][j] = m[j][i]
			else:
				m[i][j] = float(dist(d[i], d[j])) / mlen

	im = plt.imshow(m)
	im.set_cmap('autumn')
	plt.colorbar()
	plt.savefig(sys.argv[1] + ".png", bbox_inches='tight')
	plt.close()
