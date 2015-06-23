#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog <REF> <BAM> <INT q> <INT n> 

Discover and print (stdout) Context-Specific Sequencing Errors.

<REF>        Reference genome 
<BAM>        BAM file
<INT q>      Length of <q>-gram to search for
<INT n>      Maximal allowed number of Ns that the <q>-gram contains

For more details, see:

Manuel Allhoff, Alexander Schoenhuth, Marcel Martin, Ivan G. Costa, 
Sven Rahmann and Tobias Marschall.
Discovering motifs that induce sequencing errors. 
BMC Bioinformatics, 2013. 14(Suppl 5):S1, 
DOI: 10.1186/1471-2105-14-S5-S1, 
http://www.biomedcentral.com/1471-2105/14/S5/S1

Author: Manuel Allhoff
"""

from __future__ import print_function
from optparse import OptionParser
import math, sys, HTSeq, pysam, re
import scipy.misc as sc
import rpy2.robjects as robjects
import array
import pickle
import time
import copy
import scipy.stats as sps
import numpy as np
DEBUG_INDEL         = False
DEBUG_FISHER        = False
DEBUG_ASSERTIONS    = False


class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def get_annotate_qgram(genome, genome_annotate, q, search_pos):
    """Compute for each q-gram in the genome its (composed) strand bias table. 
    Consider therefore the q-gram as well as its reverse complement."""
    """#barni returns a dictionary with a 2x2 table for each qgram"""
    #consider separately pileups of q-grams last and first positions
    qgram_counts    = {}
    qgram_first     = {}
    qgram_ins_first = {}    #insertions
    qgram_del_first = {}    #deletions
    qgram_annotate      = {}
    qgram_ins_annotate  = {}
    qgram_del_annotate  = {}
    j = 0 #counter for status info
    k = 0
    l = 0
    len_genome  = len(genome)
    #pass through entire genome to analyse each q-grams' pileup
    #don't start at 0 because we're only looking for positions #search_pos before / after motif
    for i in range(search_pos, len_genome - q + 1 - search_pos):
        j += 1
        if j % 20000000 == 0: 
            print('%s / %s positions considered for q-gram annotation' %(j, len(genome)), file=sys.stderr)
        
        qgram = genome [i : i+q]
        qgram = qgram.upper()
        
        if len(qgram)!=qgram.count('A') + qgram.count('C') + qgram.count('G') + qgram.count('T'):
            #print("Warning: q-gram contains other letters than A,C,G and T, ignore q-gram",
            #        file=sys.stderr)
            #print(qgram, file=sys.stderr)
            k += 1
            continue
        else:
            l += 1

        qgram_counts[qgram] = qgram_counts[qgram] + 1 if qgram_counts.has_key(qgram) else 1
        #barni: #fm, #fmm, ... on the last position of the qgram
        #barni: list of 4 elem = 2x2 table for position i+q
        offset     = 1 #proper indexing on forward strands
        genome_pos = i + q - offset + search_pos #base position in genome to build cont. table
        qgram_effect_last = [
                genome_annotate[0][genome_pos], 
                genome_annotate[1][genome_pos],
                genome_annotate[2][genome_pos], 
                genome_annotate[3][genome_pos] ]
        #indels forward direction
        qgram_ins_effect_last = [
                genome_annotate[4][genome_pos], 
                genome_annotate[5][genome_pos], 
                genome_annotate[6][genome_pos], 
                genome_annotate[7][genome_pos] ]
        qgram_del_effect_last = [
                genome_annotate[8] [genome_pos], 
                genome_annotate[9] [genome_pos], 
                genome_annotate[10][genome_pos], 
                genome_annotate[11][genome_pos] ]

        #q-grams on forward direction, analyse therefore their last positions
        if qgram_annotate.has_key(qgram):
            qgram_annotate[qgram]     = _add_listelements(qgram_annotate[qgram], qgram_effect_last)  
            qgram_ins_annotate[qgram] = _add_listelements(qgram_ins_annotate[qgram],qgram_ins_effect_last)
            qgram_del_annotate[qgram] = _add_listelements(qgram_del_annotate[qgram],qgram_del_effect_last) 
        else:
            qgram_annotate[qgram]     = qgram_effect_last
            qgram_ins_annotate[qgram] = qgram_ins_effect_last
            qgram_del_annotate[qgram] = qgram_del_effect_last
        
        #q-gram on reverse strand, analyse therefore their first positions. 
        #Furthermore, switch read direction
        genome_pos = i - search_pos #base position in genome to build cont. table
        qgram_effect_first = [
                genome_annotate[1][genome_pos], 
                genome_annotate[0][genome_pos], 
                genome_annotate[3][genome_pos], 
                genome_annotate[2][genome_pos] ]
        #indels reverse direction
        qgram_ins_effect_first = [
                genome_annotate[5][genome_pos], 
                genome_annotate[4][genome_pos], 
                genome_annotate[7][genome_pos], 
                genome_annotate[6][genome_pos] ]
        qgram_del_effect_first = [
                genome_annotate[9] [genome_pos], 
                genome_annotate[8] [genome_pos], 
                genome_annotate[11][genome_pos], 
                genome_annotate[10][genome_pos] ]

        if qgram_first.has_key(qgram):
            qgram_first[qgram] = _add_listelements(qgram_first[qgram], qgram_effect_first)  
            qgram_ins_first[qgram] = _add_listelements(qgram_ins_first[qgram], 
                                                       qgram_ins_effect_first) 
            qgram_del_first[qgram] = _add_listelements(qgram_del_first[qgram], 
                                                       qgram_del_effect_first) 
        else:
            qgram_first[qgram] = qgram_effect_first
            qgram_ins_first[qgram] = qgram_ins_effect_first
            qgram_del_first[qgram] = qgram_del_effect_first

    #combine the q-grams on the forward and reverse strand

    # indel keys are subset of the other (normal) qgrams, only problems 
    # at beginning and end of genome
    for qgram_rev in qgram_first:
        qgram = reverse_complement(qgram_rev)

        if qgram in qgram_annotate:
            result = qgram_annotate[qgram][:]
            result = _add_listelements(result, qgram_first[qgram_rev])
            qgram_annotate[qgram] = result

            try: #insertions
                result = qgram_ins_annotate[qgram][:]
                result = _add_listelements(result, qgram_ins_first[qgram_rev])
                qgram_ins_annotate[qgram] = result
            except KeyError:
                print("Key error for reverse qgram ", qgram_rev, "for insert annotation", file=sys.stderr)
            try: #deletions
                result = qgram_del_annotate[qgram][:]
                result = _add_listelements(result, qgram_del_first[qgram_rev])
                qgram_del_annotate[qgram] = result
            except KeyError:
                print("Key error for reverse qgram ", qgram_rev, "for deletion annotation", file=sys.stderr)
        else:
            qgram_annotate[qgram]     = qgram_first[qgram_rev]
            qgram_ins_annotate[qgram] = qgram_ins_first[qgram_rev]
            qgram_del_annotate[qgram] = qgram_del_first[qgram_rev]

    print("Warning: %s q-grams of %s contain other letters than A,C,G and T, ignore these \
            q-grams" % (k, l) ,file=sys.stderr)
    return (qgram_annotate, qgram_ins_annotate, qgram_del_annotate, qgram_counts)


def reverse_complement(s, rev=True):
    """Return the reverse complement of a DNA sequence s"""
    _complement = dict(A="T", T="A", C="G", G="C", N="N")
    t = reversed(s) if rev else s
    try:
        rc = (_complement[x] for x in t)  # lazy generator expression
    except KeyError:
        return s
        
    return "".join(rc)  # join the elements into a string


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


def get_genome(ref_path, learn_chrom):
    """Return genome for random access"""
    seq = None
    i = 0
    for s in HTSeq.FastaReader(ref_path):
        if s.name != learn_chrom:
            i += 1
            continue
        else:
            seq = str(s)
            break

    if i == 1:
        learn_chrom = s.name
        seq = str(s)
    elif seq is None:
        parser.error("Sorry, the Chromosome that is using for training (%s) is not contained \
        in the reference genome (%s)! Please use -c option!" %(learn_chrom, ref_path))
    
    return seq, learn_chrom


def get_annotate_genome(genome, bampath, learn_chrom):
    """Return strand bias table for each genome position on the base of the alignment. 
    The table is represented as a list: 
    [forward match, reverse match, forward mismatch, reverse mismatch]."""
    cigar_codes = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6}
    samfile = pysam.Samfile(bampath, "rb")

    #initialize
    fm = array.array('H',[0])
    fmm = array.array('H',[0])
    rm = array.array('H',[0])
    rmm = array.array('H',[0])
    #insertions
    f_ins_m = array.array('h',[0])  # insertion to the ref
    f_ins_mm = array.array('h',[0]) # no insertion to the ref
    r_ins_m = array.array('h',[0])
    r_ins_mm = array.array('h',[0])
    #deletions
    f_del_m = array.array('H',[0])
    f_del_mm = array.array('H',[0])
    r_del_m = array.array('H',[0])
    r_del_mm = array.array('H',[0])
    
    len_genome = len(genome)
    fm          *= len_genome 
    fmm         *= len_genome 
    rm          *= len_genome 
    rmm         *= len_genome 
    f_ins_m     *= len_genome 
    f_ins_mm    *= len_genome
    r_ins_m     *= len_genome
    r_ins_mm    *= len_genome
    f_del_m     *= len_genome
    f_del_mm    *= len_genome
    r_del_m     *= len_genome
    r_del_mm    *= len_genome
    
    j = 0 #counter for status info
    #consider each read
    for read in samfile.fetch():
        if read.is_unmapped or samfile.getrname(read.rname) != learn_chrom:
            continue
    #print(samfile.getrname(read.rname), file=sys.stderr)
        j += 1
        ref_pos = read.pos
        if j % 1000000 == 0: 
            print('%s reads considered for genome annotation ' %j, file=sys.stderr)
        
        #analyse CIGAR string to consecutively compute strand bias tables
        if read.cigar is not None:
            bias, current_pos_ref, current_pos_read = 0, 0, 0 
            for code, length in read.cigar:
                if code is cigar_codes['S']:
                    current_pos_read += length
                elif code is cigar_codes['M']:
                    if read.is_reverse:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            if read.seq[i + current_pos_read + bias] == genome[pos]:
                                rm[pos]     += 1
                            else:
                                rmm[pos]    += 1
                            r_del_mm[pos]   += 1 #mismatch deletion
                            #mismatch insertion, corrected later if necessary
                            r_ins_mm[pos]   += 1 
                    else:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            if read.seq[i + current_pos_read + bias] == genome[pos]:
                                fm[pos]     += 1
                            else:
                                fmm[pos]    += 1
                            f_del_mm[pos]   += 1
                            f_ins_mm[pos]   += 1
                    bias += length
                elif code is cigar_codes['I']: 
                    # we only care about the first indel, so no loop required
                    # pos = last position in the genome where the read had M or D
                    if read.is_reverse:
                        pos = ref_pos + bias + current_pos_ref 
                        r_ins_m[pos]    += 1
                        # insert actually belongs to previous pos so delete error
                        try:
                            #only if insertion is not the first
                            if bias != 0 or current_pos_ref != 0: 
                                r_ins_mm[pos+1] -= 1 
                        except OverflowError:
                            print("Overflow. rev DEtails below", file=sys.stderr)
                            print(str(pos),str(r_ins_mm[pos-1]), str(r_ins_mm[pos+1]), str(r_ins_mm[pos]), file=sys.stderr)
                            print(read.cigar, file=sys.stderr)
                            exit(0)
                    else:
                        pos = ref_pos + bias + current_pos_ref
                        f_ins_m[pos]    += 1
                        if pos == 100391:
                            print("Cigar fwd I: ", read.cigar, file=sys.stderr) 
                        try:
                            if bias != 0 or current_pos_ref != 0: 
                                f_ins_mm[pos+1] -= 1
                        except OverflowError:
                            print("Overflow. DEtails below", file=sys.stderr)
                            print(str(pos), str(f_ins_mm[pos]), str(f_ins_mm[pos]), file=sys.stderr)
                            print(read.cigar, file=sys.stderr)
                            exit(0)
                    current_pos_read += length #manuel

                elif code is cigar_codes['D']: 
                    if read.is_reverse:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            r_del_m[pos]    += 1
                            r_ins_mm[pos]   += 1 #still mismatch for
                    else:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            f_del_m[pos]    += 1
                            f_ins_mm[pos]   += 1
                    current_pos_ref += length #manuel
                else:
                    print(code, length, file=sys.stderr)
    
    if DEBUG_ASSERTIONS:
        print("Starting debug assertions...", file=sys.stderr)
        for i in range(len_genome):
            if (f_ins_m[i] <= 0 and f_ins_mm[i] <= 0) or (r_ins_m[i] <= 0 and r_ins_mm[i] <= 0):
                print("Assertion INSERTION failed: match = mismatch = 0. Genome pos:",str(i),
                        "(f_ins_m, f_ins_mm, r_ins_m, r_ins_mm): ",
                        str(f_ins_m[i]),str(f_ins_mm[i]),str(r_ins_m[i]),str(r_ins_mm[i]), file=sys.stderr)

            if (f_del_m[i] <= 0 and f_del_mm[i] <= 0) or (r_del_m[i] <= 0 and r_del_mm[i] <= 0):
                print("Assertion DELETION failed: match = mismatch = 0. Genome pos:",str(i),
                        "(f_del_m, f_del_mm, r_del_m, r_del_mm): ",
                        str(f_del_m[i]),str(f_del_mm[i]),str(r_del_m[i]),str(r_del_mm[i]), file=sys.stderr)

            if (fm[i] == 0 and fmm[i] == 0) or (rm[i] == 0 and rmm[i] == 0):
                print("Assertion SNP (original fm, fmm, rm, rmm) failed: match = mismatch = 0. Genome pos:",str(i),
                        "(fm, fmm, rm, rmm): ", str(fm[i]),str(fmm[i]),str(rm[i]),str(rmm[i]), file=sys.stderr)
        #print("All debug assertions passed!", file=sys.stderr)
        #exit(-1)

    return (fm, rm, fmm, rmm, 
            f_ins_m, r_ins_m, f_ins_mm, r_ins_mm, 
            f_del_m, r_del_m, f_del_mm, r_del_mm)


def add_n(qgram_annotate, qgram_ins_annotate, qgram_del_annotate, n, q):
    """Extend qgram_annotate (+ qgrams for indels) by adding q-grams which contain Ns."""
    to_add      = {}
    to_add_ins  = {}
    to_add_del  = {}
    
    if n == 0:
        print('No q-grams with Ns to add' , file = sys.stderr)
        return (qgram_annotate, qgram_ins_annotate, qgram_del_annotate)

    i = 0 #counter for status info    
    #consider each possible q-gram for the given length q and number n
    for qgram_with_n in _get_all_qgrams(['A','C','G','T','N'], ['A','C','G','T','N'], q, 1, n):
        qgram_with_n = qgram_with_n[0] #q-gram that contain at least one and at most n Ns
        i += 1
        if i % 20000000 == 0: 
            print(' %s q-grams with N considered' %(i), file = sys.stderr)

        #compute all concrete q-grams of n_qgram
        possible_qgrams = get_qgramlist(qgram_with_n)
        sb_table     = [0,0,0,0] #initialize strand bias table
        sb_table_ins = [0,0,0,0] 
        sb_table_del = [0,0,0,0]
        for p_qgram in possible_qgrams:
            #each qgram key is present in all three annotations (snp, ins, del)
            if qgram_annotate.has_key(p_qgram):
                sb_table     = _add_listelements(sb_table, qgram_annotate[p_qgram][:])
                sb_table_ins = _add_listelements(sb_table_ins, qgram_ins_annotate[p_qgram][:])
                sb_table_del = _add_listelements(sb_table_del, qgram_del_annotate[p_qgram][:])
        
        #does n containing q-gram corresponds to a combosed strand bias table?
        if sb_table != [0,0,0,0]: 
            to_add[qgram_with_n] = sb_table
        if sb_table_ins != [0,0,0,0]: 
            to_add_ins[qgram_with_n] = sb_table_ins
        if sb_table_del != [0,0,0,0]: 
            to_add_del[qgram_with_n] = sb_table_del

    #extend qgram_annotate
    qgram_annotate.update(to_add)
    qgram_ins_annotate.update(to_add_ins)
    qgram_del_annotate.update(to_add_del)
    
    return (qgram_annotate, qgram_ins_annotate, qgram_del_annotate)


def _add_listelements(a, b):
    """Add i-th element of list a and list b"""
    for i in range(len(a)):
        a[i] += b[i]
    
    return a


def get_qgramlist(qgram_with_n):
    """Return list of all possible q-grams for q-grams which contain Ns"""
    if qgram_with_n.count("N")==0:
        return [qgram_with_n]
    return _get_qgramlist_help(qgram_with_n, [])


def _get_qgramlist_help(qgram, result):
    """Substitute N with A, C, G or T and return resulting q-grams"""
    if qgram.count('N') == 0:
        result.append(qgram)
    else:
        i = qgram.index("N")
        for n in ['A','C','G','T']:
            new_qgram = qgram[:i] + n + qgram[i+1:]
            _get_qgramlist_help(new_qgram, result)
        return result


def _get_all_qgrams(alphabet, erg, length, level, n):
    """Return all possible q-grams of the given length, over the given alphabet
    and with at least one and at most n Ns"""
    if length == level:
        yield erg
    else:
        for letter in alphabet:
            for el in erg:
                #not too many Ns
                if letter == 'N' and el.count('N') >= n:
                    continue
                #not too less Ns
                if length - level <= 1 and el.count('N') == 0 and letter != 'N':
                    continue

                for r in _get_all_qgrams(alphabet, [el + letter], length, level+1, n):
                    yield r


def get_sb_score(qgram_annotate):
    """Calculate Strand Bias score (based on p-value) for each q-gram"""
    results = []
    i = 0 #counter for status info
    print('Start Strand Bias Score calculation', file=sys.stderr)
    for k in qgram_annotate.keys():
        i += 1
        if i % 1000000 == 0: 
            print(" %s / %s Strand Bias Scores calculated" %(i, len(qgram_annotate.keys())), file=sys.stderr)
        
        #get p-value for the strand bias table of q-gram k
        p_value = get_pvalue(qgram_annotate[k][0], qgram_annotate[k][1], qgram_annotate[k][2], qgram_annotate[k][3])
        
        #compute negative logarithm (base 10) of p-value or, if necessary, set to maxint
        strand_bias_score = sys.maxint if p_value < 1/10.0**300 else -math.log10(p_value)
        
        results.append((k, qgram_annotate[k][0], qgram_annotate[k][1], qgram_annotate[k][2], qgram_annotate[k][3], strand_bias_score))
    
    return results


def output(results, genome, task, qgram_counts):
    """Output the results"""
    print("###############################",task,"#########################")
    print("#Sequence", "Occurrence", "Forward Match", "Backward Match", 
            "Forward Mismatch", "Backward Mismatch", "Strand Bias Score", "FER (Forward Error Rate)", 
            "RER (Reverse Error Rate), ERD (Error rate Difference)", sep = '\t')
    for seq, f_m, r_m, f_mm, r_mm, sb_score, fer, rer, erd in results:
        occ = count_app(seq, qgram_counts)
        print(seq, occ, f_m, r_m, f_mm, r_mm, sb_score, fer, rer, erd, sep = '\t')


def get_motifspace_size(q,n):
    """return length of search space according to equation which is mentioned in Section 3.1 of the paper"""
    return reduce(lambda x, y: x + (int(sc.comb(q, y, exact=True)) * 4**(q-y)), 
            [i for i in range(1, n+1)], int(sc.comb(q,0,exact=True)) * 4**(q-0))

#new count function
def count_app(qgram_, qgram_counts):
    res = 0
    rev = reverse_complement(qgram_)
    #forward
    for qgram in get_qgramlist(qgram_):
        if qgram_counts.has_key(qgram):
            res += qgram_counts[qgram]
    #backward
    for qgram in get_qgramlist(rev):
        if qgram_counts.has_key(qgram):
            res += qgram_counts[qgram]

    return res


def count(qgram, genome):
    """Count number of q-grams and its reverse complement in genome"""
    rev = reverse_complement(qgram)
    rev = rev.replace('N', '.')
    qgram = qgram.replace('N', '.')
    
    return  len([m.start() for m in re.finditer(r'(?=(%s))' %qgram, genome)] + [m.start() for m in re.finditer(r'(?=(%s))' %rev, genome)])


def log(results, s, genome):
    """Result is a list [(seq, fm, rm, fmm, rmm, p_value)] """
    if file_indel_log == "":
        return
    try:
        with open(file_indel_log, 'a') as f:
            print("#Sequence", "Occurrence", "Forward Match", "Backward Match", 
                  "Forward Mismatch", "Backward Mismatch", "Strand Bias Score",
                  "FER","RER", sep = '\t', file = f)
            for r in results:
                seq = r[0]
                fer = float(r[3]) / (r[1] + r[3]) #forward error rate
                rer = float(r[4]) / (r[2] + r[4]) #forward error rate
                print(seq,count(seq,genome),str(r[1]),str(r[2]),str(r[3]),
                        str(r[4]),str(r[5]),str(fer),str(rer),sep='\t',file=f)
    except IOError:
        print("Could not open indel log file.", file=sys.stderr)


def ident(genome, genome_annotate, q, n, alpha=0.05, epsilon=0.03, delta=0.05, search_pos=0, only_indels=False):
    """Identify critical <q>-grams (with <n> Ns) with reference to significance and error rate""" 
    motifspacesize_log = math.log10(get_motifspace_size(q, n))
    alpha_log = math.log10(float(alpha))

    #annotate each q-gram with Strand Bias Table
    qgram_annotate, qgram_ins_annotate, qgram_del_annotate, qgram_counts = \
            get_annotate_qgram(genome, genome_annotate, q, search_pos) 
    print("Number of in-, del keys before all_results [should match!]: ", len(qgram_annotate), 
            len(qgram_ins_annotate), len(qgram_del_annotate), file=sys.stderr)

    #extend set of q-grams with q-grams containing Ns
    add_n(qgram_annotate, qgram_ins_annotate, qgram_del_annotate, n, q) 
    #annotate each q-gram with Strand Bias Score
    all_results = (get_sb_score(qgram_annotate),
            get_sb_score(qgram_ins_annotate),
            get_sb_score(qgram_del_annotate))
    print("Number of in-, del keys after all_results [should match!]: ",
            len(all_results[1]),len(all_results[1]), len(all_results[2]), file=sys.stderr)
    #log indels
    log(all_results[1], "insertions", genome)
    log(all_results[2], "deletions", genome)
    #filter statistically significant motifs (Bonferroni Correction)
    all_sig_results = (filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results[0]),
            filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results[1]), 
            filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results[2])) 

    tasks = { 0 : 'snps', 1 : 'insertion', 2 : 'deletions' }
    for index, sig_results in enumerate(all_sig_results):
        # skip SNPs output
        # TODO: this should be done in a more generic way, so as to avoid calculating
        # SNP-inducing motifs at all
        if index == 0 and only_indels:
            continue

        #filter motifs with regards to background error rate (epsilon) error rate difference (delta)
        results = []
        for seq, f_m, r_m, f_mm, r_mm, p_value_score in sig_results:
            try:
                fer = float(f_mm) / (f_m + f_mm) #forward error rate
            except ZeroDivisionError:
                print("ZeroDivisionError fer:",count_app(seq,genome),tasks[index],
                        seq,str(f_m),str(r_m),str(f_mm),str(r_mm),str(p_value_score), file=sys.stderr)
                continue
            try:
                rer = float(r_mm) / (r_m + r_mm) #reverse error rate
            except ZeroDivisionError:
                print("ZeroDivisionError rer:",count(seq,genome),tasks[index],
                        seq,str(f_m),str(r_m),str(f_mm),str(r_mm),str(p_value_score), file=sys.stderr)
                continue

            #TODO remove, debug output
            if index > 0 and DEBUG_INDEL is True:
                print("(",tasks[index],count(seq,genome),seq,str(f_m),str(r_m),
                        str(f_mm),str(r_mm),str(p_value_score),str(fer),str(rer),")", file=sys.stderr)

            if rer < epsilon: #filter sequences with too high epsilon (background error rate)
                erd = fer - rer #error rate difference
                if erd >= delta: #filter sequences with too low delta (error rate difference cutoff)
                    results.append((seq, f_m, r_m, f_mm, r_mm, p_value_score, fer, rer, erd)) 
            
        results.sort(key=lambda x: x[8],reverse=True) #sort by erd (error rate difference)
        output(results, genome, tasks[index], qgram_counts)


if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    
    parser.add_option("-a", dest="alpha", default=0.05, type="float", help="FWER (family-wise error rate) alpha, default: 0.05")
    parser.add_option("-e", dest="epsilon", default=0.03, type="float", help="background error rate cutoff epsilon, default: 0.03")
    parser.add_option("-d", dest="delta", default=0.05, type="float", help="error rate difference cutoff delta, default: 0.05")
    parser.add_option("-c", dest="learn_chrom", default="chr1", help="chromosome that is used to derive Context Specific Errors, default: chr1")
    parser.add_option("-v", dest="version", default=False, action="store_true", help="show script's version")
    parser.add_option("-p", dest="position", default=0, help="search for motif effect at position p after the motif. Default is 0 (last position of motif)")
    parser.add_option("--log-indels", dest="log_indels", default="", help="log indels in file")
    parser.add_option("--verbose", dest="verbosity_", default=0, help="print verbose messages")
    parser.add_option("-s", dest="serialize", default="", help="dump genome annotation, ie contingency tables for genome positions")
    parser.add_option("-l", dest="load", default="", help="load serialized object instead of creating contingency tables for genome positions")
    parser.add_option("--only-indels", dest="only_indels", default=False, action="store_true", help="only output indels-inducing motifs, not SNPs")
    
    (options, args) = parser.parse_args()
    
    if options.version:
            version = "version \"0.32\""
            print("")
            print(version)
            sys.exit()
    
    if len(args) != 4:
        parser.error("Sorry, exactly four parameters are required.")
        parser.print_help(sys.stderr)
    
    #map arguments
    refpath = args[0]
    bampath = args[1]
    q = int(args[2])
    n = int(args[3])

    genome, options.learn_chrom = get_genome(refpath, options.learn_chrom)

    if options.load != "":
        #load input from serialized file
        print("Loading serialized object instead of annotating (parsing) genome...",file=sys.stderr)
        with open(options.load, 'rb') as inpkl:
            genome_annotate = pickle.load(inpkl)
    else:
        print("Annotating (parsing) genome...", file=sys.stderr)
        genome_annotate = get_annotate_genome(genome, bampath, options.learn_chrom)

    if options.serialize != "":
        #object serialization
        print("Dumping serialized object in 3 seconds...", file=sys.stderr)
        time.sleep(3)
        with open(options.serialize, 'wb') as outpkl:
            pickle.dump(genome_annotate, outpkl, pickle.HIGHEST_PROTOCOL)

    global file_indel_log, verbosity
    verbosity       = options.verbosity_
    file_indel_log  = options.log_indels

    if file_indel_log != "":
        print("Warning! Dumping indel details will take a long time...",file=sys.stderr)

    print("Searching position is motif[last + %s]" % options.position, file=sys.stderr) 
    print("#Searching position is motif[last + %s]" % options.position)
    ident(genome, genome_annotate, q, n, options.alpha, options.epsilon, 
            options.delta, int(options.position),options.only_indels)
