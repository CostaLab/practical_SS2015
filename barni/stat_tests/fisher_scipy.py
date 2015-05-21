#!/usr/bin/python

from __future__ import print_function
from optparse import OptionParser
import math, sys, HTSeq, pysam, re
import rpy2.robjects as robjects
from fisher import pvalue
import array
import pickle
import copy
import scipy.stats as sps
import numpy as np

forward_match   = 2874
reverse_match   = 3871
forward_mismatch= 157
reverse_mismatch= 2

#scipy chi squared test
arr         = np.array([[forward_match, forward_mismatch],[reverse_match, reverse_mismatch]])
#scipy fisher's exact test
for i in range (10000):
    sc_odds, sc_fish_pv = sps.fisher_exact(arr)

s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n\nScipy Fisher's Exact\t(2-sided)"
s += "\n"
s += "\t\t".join(["","", str(sc_fish_pv)])
s += "\n==================="
print(s)
