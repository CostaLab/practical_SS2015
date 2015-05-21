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

forward_match   = 15493
reverse_match   = 17481
forward_mismatch= 15747
reverse_mismatch= 26000

ff      = robjects.r['fisher.test']
fchi    = robjects.r['chisq.test']

matrix  = [forward_match, reverse_match, forward_mismatch, reverse_mismatch]
table   = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)
table2  = copy.deepcopy(table)

# R values
rf_pvalue    = tuple(ff(table)[0])[0]
rchi_pvalue  = tuple(fchi(table2)[2])[0]

# new fisher's test library
libf_pvalue1 = pvalue(forward_match, forward_mismatch, reverse_match, reverse_mismatch)
libf_pvalue2 = pvalue(forward_match, reverse_match, forward_mismatch, reverse_mismatch)

#scipy chi squared test
arr         = np.array([[forward_match, forward_mismatch],[reverse_match, reverse_mismatch]])
sc_chi2, sc_p, sc_dof, sc_expected = sps.chi2_contingency(arr)

#scipy fisher's exact test
sc_odds_l, sc_fish_pv_left = sps.fisher_exact(arr, alternative="less")
sc_odds_r, sc_fish_pv_right = sps.fisher_exact(arr, alternative="greater")
sc_odds, sc_fish_pv = sps.fisher_exact(arr, alternative="two-sided")

s = "-------------------\n" 
s += "fm, rm, fmm, rmm = " + ", ".join([str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch)])
s += "\n"
s += "\t\t\t".join(["R test results:","fisher", "chi"])
s += "\n"
s += "\t\t ".join(["","",str(rf_pvalue), str(rchi_pvalue)])
s += "\n\n"
s += "\t\t".join(["New library fisher's:","fisher.py [left]" + "\tfisher.py [right]" + "\tfisher.py [two]"])
s += "\n"
s += "\t\t ".join(["","",str(libf_pvalue1.left_tail), str(libf_pvalue1.right_tail), str(libf_pvalue1.two_tail)])
s += "\n\nNew library fisher's with (fm, rm, fmm, rmm) instead of (fm, fmm, rm, rmm), like above \n" 
s += "\t\t ".join(["","",str(libf_pvalue2.left_tail), str(libf_pvalue2.right_tail), str(libf_pvalue2.two_tail)])
s += "\n\nScipy chi-squared\t\t" + "\t\t".join(["chi2", "\tp_value", "\tdof (degrees of freedom)"])
s += "\n"
s += "\t\t ".join(["","",str(sc_chi2), str(sc_p), str(sc_dof)])
s += "\n\nScipy Fisher's Exact\t\tp-value left\tright\t(2-sided)"
s += "\n"
s += "\t\t\t\t".join(["","",str(sc_fish_pv_left), str(sc_fish_pv_right), str(sc_fish_pv)])
s += "\n==================="
print(s)
