#!/usr/bin/python

from __future__ import print_function
import pysam
import sys
import os
import scipy.misc as sc
from collections import defaultdict

def get_motifspace_size(q,n):
    """return length of search space according to equation which is mentioned in Section 3.1 of the paper"""
    return reduce(lambda x, y: x + (int(sc.comb(q, y, exact=True)) * 4**(q-y)), 
            [i for i in range(1, n+1)], int(sc.comb(q,0,exact=True)) * 4**(q-0))

if len(sys.argv) < 3:
	print("ERROR:", os.path.basename(sys.argv[0]), "FASTA k", file = sys.stderr)
	exit(1)

if not os.path.isfile(sys.argv[1]):
	print("ERROR: FASTA must be a valid file", file = sys.stderr)
	exit(1) 

k = int(sys.argv[2])

if k < 1:
	print("ERROR: k must be >= 1", file = sys.stderr)
	exit(1)

def reverse_complement(s, rev=True):
    """Return the reverse complement of a DNA sequence s"""
    _complement = dict(A="T", T="A", C="G", G="C", N="N")
    t = reversed(s) if rev else s

    rc = (_complement[x] for x in t)  # lazy generator expression
        
    return "".join(rc)  # join the elements into a string

ffile = pysam.Fastafile(sys.argv[1])

s = set()

max_kmers = get_motifspace_size(k, 0)

for ref in ffile.references:
	seq = ffile.fetch(reference=ref)

	length = len(seq)

	print("# calculating %d-mers for %s.." % (k, ref))
	print("# seq length:", length)
	revseq = reverse_complement(seq)
	
	for i in range(len(seq) - k + 1):
		motif = seq[i:i+k]

		if "N" in motif:
			continue

		s.add(motif)

	for i in range(len(revseq) - k + 1):
		motif = revseq[i:i+k]

		if "N" in motif:
			continue

		s.add(motif)

ffile.close()

tot_kmers = len(s)
perc = (float(tot_kmers) / float(max_kmers)) * 100

print("# found %d-mers: %d out of %d (%f%%)" % (k, tot_kmers, max_kmers, perc))
print("# dumping to file in lexicographic order")

for motif in sorted(s):
	print(motif)
