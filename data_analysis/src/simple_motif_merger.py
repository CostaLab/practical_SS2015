#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
%prog [options] <INT q> <INT n> <FILE> <FILE> ... <FILE>

Takes a list of .data files (from the CSE discovery script) and,
for common motifs, recalculates the congintency tables, applies
the cutoff thresholds and outputs the new motif data.

You must provide at least two files.

<INT q>      Length of <q>-gram to search for
<INT n>      Maximal allowed number of Ns that the <q>-gram contains
<FILE>       Tab-separated data file from discovery_cse script

Author: Fabio Ticconi
Based on code from: Manuel Allhoff
"""

from __future__ import print_function
from optparse import OptionParser

import sys
import itertools
import scipy.stats as sps
import scipy.misc as sc
import numpy as np
import rpy2.robjects as robjects
import math

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def load_file(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            if line[0] == '#' or line[0] == ' ':
                continue

            # Motif, Occurrences, FM, BM, FMM, BMM, SBS, FER, RER, ERD
            fields = line.split()

            motif = fields[0]
            opts  = [int(x) for x in fields[1:6]]

            # accumulate tables of same motif (there shouldn't be any here,
            # unless the file was a naive merge of multiple results files)
            if motif in d:
                d[motif].append(opts)
            else:
                d[motif] = [opts]

    # if there are any duplicate motifs, we merge their tables together
    for l in d:
        d[l] = [sum(sb) for sb in itertools.izip(*(d[l]))]

    return d

def merge_dicts(list_of_dicts, d):
    if not list_of_dicts:
        return d
    
    d2 = list_of_dicts[0]
    
    d3 = {k:[v, d2[k]] for k,v in d.iteritems() if k in d2}
    d = {k:[sum(sb) for sb in itertools.izip(*(v))] for k,v in d3.iteritems()}

    return merge_dicts(list_of_dicts[1:], d)
    
def get_motifspace_size(q,n):
    """return length of search space according to equation which is mentioned in Section 3.1 of the paper"""
    return reduce(lambda x, y: x + (int(sc.comb(q, y, exact=True)) * 4**(q-y)), 
            [i for i in range(1, n+1)], int(sc.comb(q,0,exact=True)) * 4**(q-0))
    
def get_pvalue(fm, rm, fmm, rmm):
    """Return p-value of given Strand Bias Table"""
    #decide whether Fisher's exact Test or ChiSq-Test should be used
    limit = 5000
    if (fm > limit or rm > limit or fmm > limit or rmm > limit):
        #for correct chi squared test
        if ((fm > 0 or fmm > 0) and (rm > 0 or rmm > 0) 
                and (fm > 0 or rm > 0) and (fmm > 0 or rmm > 0)): 
            arr = np.array([[fm, fmm],[rm, rmm]])
            try:
                sc_chi2, p_value, sc_dof, sc_expected = sps.chi2_contingency(arr)
                return p_value
            except ValueError:
                print("Chi calculation error. fm, rm, fmm, rmm: ",str(fm), str(rm),
                        str(fmm),str(rmm), file=sys.stderr)
                exit(-1)
        else:
            return 1.0

    f = robjects.r['fisher.test']
    matrix = [fm, rm, fmm, rmm]
    table = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)
    p_value_tmp = f(table)[0] 
    p_value = tuple(p_value_tmp)[0] #some necessary magic for r object

    return p_value

def output(results):
    """Output the results"""
    print("#Sequence", "Occurrence", "Forward Match", "Backward Match", "Forward Mismatch",
          "Backward Mismatch", "Strand Bias Score", "FER (Forward Error Rate)",
          "RER (Reverse Error Rate), ERD (Error rate Difference)", sep = '\t')
    
    for seq, occ, forward_match, reverse_match, forward_mismatch, reverse_mismatch, sb_score, fer, rer, erd in results:
        print(seq, occ, forward_match, reverse_match, forward_mismatch, reverse_mismatch, sb_score, fer, rer, erd, sep = '\t')

if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    
    parser.add_option("-a", dest="alpha", default=0.05, type="float", help="FWER (family-wise error rate) alpha, default: 0.05")
    parser.add_option("-e", dest="epsilon", default=0.03, type="float", help="background error rate cutoff epsilon, default: 0.03")
    parser.add_option("-d", dest="delta", default=0.05, type="float", help="error rate difference cutoff delta, default: 0.05")
    
    (options, args) = parser.parse_args()


    if len(sys.argv) < 5:
        parser.error("must provide at least four arguments")
        parser.print_help(sys.stderr)

    q = int(sys.argv[1])
    n = int(sys.argv[2])
    files = sys.argv[3:]

    list_of_dicts = [load_file(f) for f in files]

    d = merge_dicts(list_of_dicts[1:], list_of_dicts[0])

    # FIXME: is this correct? Or are, now, the "hypothesis"
    # as many as the files (as opposed to the search space)?
    bnf = alpha / get_motifspace_size(q, n)

    #print("Bonferroni threshold:", bnf, file=sys.stderr)

    results = []
    for motif, fields in d.iteritems():
        occ, fm, rm, fmm, rmm = fields[0:5]

        pv = get_pvalue(fm, rm, fmm, rmm)
        
        #compute negative logarithm (base 10) of p-value or, if necessary, set to maxint
        sbs = sys.maxint if pv < 1/10.0**300 else -math.log10(pv)

        fer = float(fmm) / (fm + fmm) #forward error rate
        rer = float(rmm) / (rm + rmm) #reverse error rate
        erd = fer - rer #error rate difference

        # Bonferroni correction
        if pv > bnf:
            #print("Pvalue too high:", motif, occ, fm, rm, fmm, rmm, sbs, fer, rer, erd, file=sys.stderr, sep='\t')
            continue

        # background error rate correction
        if rer >= epsilon:
            #print("RER too big:", motif, occ, fm, rm, fmm, rmm, sbs, fer, rer, erd, file=sys.stderr, sep='\t')
            continue

        if erd < delta:
            #print("ERD too small:", motif, occ, fm, rm, fmm, rmm, sbs, fer, rer, erd, file=sys.stderr, sep='\t')
            continue

        results.append((motif, occ, fm, rm, fmm, rmm, sbs, fer, rer, erd)) 
            
    results.sort(key=lambda x: x[9],reverse=True) #sort by erd (error rate difference)
    output(results)