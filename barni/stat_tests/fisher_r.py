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

ff      = robjects.r['fisher.test']
fchi    = robjects.r['chisq.test']

matrix  = [forward_match, reverse_match, forward_mismatch, reverse_mismatch]
table   = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)

# R values
for i in range (10000):
    rf_pvalue    = tuple(ff(table)[0])[0]
rchi_pvalue  = -1


s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n"
s += "\t\t\t".join(["R test results:","fisher", "chi"])
s += "\n"
s += "\t\t ".join(["","",str(rf_pvalue), str(rchi_pvalue)])
s += "\n==================="
print(s)
