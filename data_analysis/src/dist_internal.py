#!/usr/bin/python

import sys
from ctypes import *

if len(sys.argv) != 2:
	print "must provide exactly 1 .data file"
	exit(1)

lib = cdll.LoadLibrary("./libdistint.so.1.0")

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

array   = (c_char_p * len(d))(*d)
array_l = c_uint(len(d))
motif_l = c_uint(mlen)

dist = lib.get_internal_distance

dist.restype = c_double

print dist(array, array_l, motif_l)
