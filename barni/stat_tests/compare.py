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

forward_match   = 13137
reverse_match   = 14372
forward_mismatch= 993
reverse_mismatch= 77

#scipy chi squared test
arr         = np.array([[forward_match, forward_mismatch],[reverse_match, reverse_mismatch]])
sc_chi2, sc_p, sc_dof, sc_expected = sps.chi2_contingency(arr)

# R
fchi    = robjects.r['chisq.test']
matrix  = [forward_match, reverse_match, forward_mismatch, reverse_mismatch]
table   = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)
# R values
rchi_pvalue  = tuple(fchi(table)[2])[0]
print("Equal? : ", str(sc_p == rchi_pvalue))

s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n"
s += "\n\nScipy chi-squared\t\t" + "\t\t".join(["chi2", "\tp_value", "\tdof (degrees of freedom)"])
s += "\n"
s += "\t\t ".join(["","",str(sc_chi2), str(sc_p), str(sc_dof)])
s += "\n==================="
print(s)

s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n"
s += "\t\t\t".join(["R test results:", "chi"])
s += "\n"
s += "\t\t ".join(["","",str(rchi_pvalue)])
s += "\n==================="
print(s)
