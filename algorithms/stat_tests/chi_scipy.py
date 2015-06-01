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

forward_match   = 0
reverse_match   = 1574
forward_mismatch= 0
reverse_mismatch= 77

#scipy chi squared test
arr         = np.array([[forward_match, forward_mismatch],[reverse_match, reverse_mismatch]])
print(arr)
for i in range (10):
    if (True) and not (forward_match == 0 and forward_mismatch == 0):
        sc_chi2, sc_p, sc_dof, sc_expected = sps.chi2_contingency(arr)
    else:
        print("shit")


s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n"
s += "\n\nScipy chi-squared\t\t" + "\t\t".join(["chi2", "\tp_value", "\tdof (degrees of freedom)"])
s += "\n"
s += "\t\t ".join(["","",str(sc_chi2), str(sc_p), str(sc_dof)])
s += "\n==================="
print(s)
