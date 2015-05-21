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

forward_match   = 154930
reverse_match   = 17481
forward_mismatch= 15747
reverse_mismatch= 260000

ff      = robjects.r['fisher.test']
fchi    = robjects.r['chisq.test']

matrix  = [forward_match, reverse_match, forward_mismatch, reverse_mismatch]
table   = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)

# R values
rf_pvalue    = -1
for i in range (100000):
    rchi_pvalue  = tuple(fchi(table)[2])[0]



s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n"
s += "\t\t\t".join(["R test results:","fisher", "chi"])
s += "\n"
s += "\t\t ".join(["","",str(rf_pvalue), str(rchi_pvalue)])
s += "\n==================="
print(s)
