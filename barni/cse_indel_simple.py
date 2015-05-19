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
from fisher import pvalue
import array
import pickle
import time
import copy
DEBUG_INDEL = False
DEBUG_FISHER= True
LOG_DEL     = "../output/bsubtilis/del_3_80.log"

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def get_annotate_qgram(genome, genome_annotate, q):
    """Compute for each q-gram in the genome its (composed) strand bias table. 
    Consider therefore the q-gram as well as its reverse complement."""
    """#barni returns a dictionary with a 2x2 table for each qgram"""
    #consider separately pileups of q-grams last and first positions
    qgram_last  = {}
    qgram_first = {}
    qgram_ins_last  = {}    #insertions
    qgram_ins_first = {}    #insertions
    qgram_del_last  = {}    #deletions
    qgram_del_first = {}    #deletions
    j = 0 #counter for status info
    k = 0
    l = 0
    len_genome  = len(genome)
    offset      = 0 #if 0 search for last position of motif,else last + offset
    #pass through entire genome to analyse each q-grams' pileup
    for i in range(len(genome) - q):
        j += 1
        if j % 20000000 == 0: 
            print('%s / %s positions considered for q-gram annotation' %(j, len(genome)), file=sys.stderr)
        
        qgram = genome [i : i+q]
        qgram = qgram.upper()
        
        if len(qgram) != qgram.count('A') + qgram.count('C') + qgram.count('G') + qgram.count('T'):
            #print("Warning: q-gram contains other letters than A,C,G and T, ignore q-gram" ,file=sys.stderr)
            #print(qgram, file=sys.stderr)
            k += 1
            continue
        else:
            l += 1

        #barni: #fm, #fmm, ... on the last position of the qgram
        #barni: list of 4 elem = 2x2 table for position i+q
        qgram_effect_last = [
                genome_annotate[0][i + q], 
                genome_annotate[1][i + q],
                genome_annotate[2][i + q], 
                genome_annotate[3][i + q] ]

        #q-grams on forward direction, analyse therefore their last positions
        if qgram_last.has_key(qgram):
            qgram_last[qgram]     = _add_listelements(qgram_last[qgram], qgram_effect_last)  
        else:
            qgram_last[qgram]     = qgram_effect_last

        ################################################
        #indels forward direction
        if i + q + offset <= len_genome:
            qgram_ins_effect_last = [
                    genome_annotate[4][i + q + offset], 
                    genome_annotate[5][i + q + offset], 
                    genome_annotate[6][i + q + offset], 
                    genome_annotate[7][i + q + offset] ]
            qgram_del_effect_last = [
                    genome_annotate[8] [i + q + offset], 
                    genome_annotate[9] [i + q + offset], 
                    genome_annotate[10][i + q + offset], 
                    genome_annotate[11][i + q + offset] ]

            if qgram_ins_last.has_key(qgram):
                qgram_ins_last[qgram] = _add_listelements(qgram_ins_last[qgram], qgram_ins_effect_last) 
                qgram_del_last[qgram] = _add_listelements(qgram_del_last[qgram], qgram_del_effect_last) 
            else:
                qgram_ins_last[qgram] = qgram_ins_effect_last
                qgram_del_last[qgram] = qgram_del_effect_last
        
        
        #q-gram on reverse strand, analyse therefore their first positions. Furthermore, switch read direction
        qgram_effect_first = [
                genome_annotate[1][i], 
                genome_annotate[0][i], 
                genome_annotate[3][i], 
                genome_annotate[2][i] ]
        if qgram_first.has_key(qgram):
            qgram_first[qgram] = _add_listelements(qgram_first[qgram], qgram_effect_first)  
        else:
            qgram_first[qgram] = qgram_effect_first

        ################################################
        #indels reverse direction
        if i - offset >= 0:
            qgram_ins_effect_first = [
                    genome_annotate[5][i - offset], 
                    genome_annotate[4][i - offset], 
                    genome_annotate[7][i - offset], 
                    genome_annotate[6][i - offset] ]
            qgram_del_effect_first = [
                    genome_annotate[9] [i - offset], 
                    genome_annotate[8] [i - offset], 
                    genome_annotate[11][i - offset], 
                    genome_annotate[10][i - offset] ]

            if qgram_ins_first.has_key(qgram):
                qgram_ins_first[qgram] = _add_listelements(qgram_ins_first[qgram], qgram_ins_effect_first) 
                qgram_del_first[qgram] = _add_listelements(qgram_del_first[qgram], qgram_del_effect_first) 
            else:
                qgram_ins_first[qgram] = qgram_ins_effect_first
                qgram_del_first[qgram] = qgram_del_effect_first
        
    #combine the q-grams on the forward and reverse strand
    qgram_annotate      = {}
    qgram_ins_annotate  = {}
    qgram_del_annotate  = {}

    # indel keys are subset of the other (normal) qgrams, only problems 
    # at beginning and end of genome
    for qgram in qgram_last: 
        qgram_annotate[qgram] = qgram_last[qgram]
        try:
            qgram_ins_annotate[qgram] = qgram_ins_last[qgram]
            qgram_del_annotate[qgram] = qgram_del_last[qgram]
        except KeyError:
            print("Key error for qgram ", qgram, "for indel annotation")
        qgram_rev = reverse_complement(qgram)
        if qgram_rev in qgram_first:
            result = qgram_annotate[qgram][:]
            result = _add_listelements(result, qgram_first[qgram_rev])
            qgram_annotate[qgram] = result

            try: #insertions
                result = qgram_ins_annotate[qgram][:]
                result = _add_listelements(result, qgram_ins_first[qgram_rev])
                qgram_ins_annotate[qgram] = result
            except KeyError:
                print("Key error for reverse qgram ", qgram, "for insert annotation")
            try: #deletions
                result = qgram_del_annotate[qgram][:]
                result = _add_listelements(result, qgram_del_first[qgram_rev])
                qgram_del_annotate[qgram] = result
            except KeyError:
                print("Key error for reverse qgram ", qgram, "for deletion annotation")

    print("Warning: %s q-grams of %s contain other letters than A,C,G and T, ignore these q-grams" %(k, l) ,file=sys.stderr)
    return (qgram_annotate, qgram_ins_annotate, qgram_del_annotate)


def reverse_complement(s, rev=True):
    """Return the reverse complement of a DNA sequence s"""
    _complement = dict(A="T", T="A", C="G", G="C", N="N")
    t = reversed(s) if rev else s
    try:
        rc = (_complement[x] for x in t)  # lazy generator expression
    except KeyError:
        return s
        
    return "".join(rc)  # join the elements into a string


def computeNewFisher(fm, fmm, rm, rmm):
    fisher_pvalue = pvalue(fm, fmm, rm, rmm) #returns (left_t, right_t, two_t)
    #print("forward",str(chi), str(fisher_pvalue.left_tail), str(fisher_pvalue.right_tail), str(fisher_pvalue.two_tail))
    return fisher_pvalue


def get_pvalue(forward_match, reverse_match, forward_mismatch, reverse_mismatch):
    """Return p-value of given Strand Bias Table"""
        #decide whether Fisher's exact Test or ChiSq-Test should be used
    limit = 5000
    if forward_match > limit or reverse_match > limit or forward_mismatch > limit or reverse_mismatch > limit:
        f = robjects.r['chisq.test']
        test = 'chisq'
    else:
        f = robjects.r['fisher.test']
        test = 'fisher'
    
    matrix = [forward_match, reverse_match, forward_mismatch, reverse_mismatch]
    table = robjects.r.matrix(robjects.IntVector(matrix), nrow=2)
    p_value_tmp = f(table)[0] if test == 'fisher' else f(table)[2]
    p_value = tuple(p_value_tmp)[0] #some necessary magic for r object

    #TODO remove
    if DEBUG_FISHER:
        table2 = copy.deepcopy(table)
        rf_pvalue    = tuple(f(table)[0])[0]
        rchi_pvalue  = tuple(f(table2)[2])[0]
        libf_pvalue1 = computeNewFisher(forward_match, forward_mismatch, reverse_match, reverse_mismatch)
        libf_pvalue2 = computeNewFisher(forward_match, reverse_match, forward_mismatch, reverse_mismatch)
        s = "-------------------\n" 
        s += "(fm, rm, fmm, rmm) (", str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch),")\n"
        s += "\t\t".join(["fisher", "chi", "fisher.py [left]", "fisher.py [right]", "fisher.py [two]"])
        s += "\n"
        s += "\t ".join([str(rf_pvalue), str(rchi_pvalue), str(libf_pvalue1.left_tail), str(libf_pvalue1.right_tail), str(libf_pvalue1.two_tail)])
        s += "\nlib fisher using fm, rm, fmm, rmm, above was (fm,fmm,rm,rmm) used\n" 
        s += "\t ".join([' ',' ', str(libf_pvalue2.left_tail), str(libf_pvalue2.right_tail), str(libf_pvalue2.two_tail)])
        s += "\n==================="
        print(s)

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
    #barni return a tuple (fm, fmm, rm, rmm) of 4 arrays
    each of these representing the genome, and for each position in the genome they contain
    the #matches and #mismatches
    """
    cigar_codes = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6}
    samfile = pysam.Samfile(bampath, "rb")

    #initialize
    fm = array.array('H')
    fmm = array.array('H')
    rm = array.array('H')
    rmm = array.array('H')
    #insertions
    f_ins_m = array.array('H')  # insertion to the ref
    f_ins_mm = array.array('H') # no insertion to the ref
    r_ins_m = array.array('H')
    r_ins_mm = array.array('H')
    #deletions
    f_del_m = array.array('H')
    f_del_mm = array.array('H')
    r_del_m = array.array('H')
    r_del_mm = array.array('H')
    
    for i in range(len(genome)):
        fm.append(0)
        fmm.append(0)
        rm.append(0)
        rmm.append(0)
        f_ins_m.append(0)
        f_ins_mm.append(0) 
        r_ins_m.append(0)
        r_ins_mm.append(0)
        f_del_m.append(0)
        f_del_mm.append(0)
        r_del_m.append(0)
        r_del_mm.append(0)
    
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
                                rm[pos] += 1
                            else:
                                rmm[pos] += 1
                            #there is a mismatch for insertion and deletion as well
                            r_del_mm[pos] += 1
                            r_ins_mm[pos] += 1
                    else:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            if read.seq[i + current_pos_read + bias] == genome[pos]:
                                fm[pos] += 1
                            else:
                                fmm[pos] += 1
                            #there is a mismatch for insertion and deletion as well
                            f_del_mm[pos] += 1
                            f_ins_mm[pos] += 1
                    bias += length
                elif code is cigar_codes['I']: 
                    if read.is_reverse:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            r_ins_m[pos] += 1
                            #still mismatch for 
                            r_del_mm[pos] += 1
                    else:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            f_ins_m[pos] += 1
                            #still mismatch for 
                            f_del_mm[pos] += 1
                    current_pos_read += length #manuel
                elif code is cigar_codes['D']: 
                    if read.is_reverse:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            r_del_m[pos] += 1
                            #still mismatch for
                            r_ins_mm[pos] += 1
                    else:
                        for i in range(length):
                            pos = ref_pos + bias + current_pos_ref + i
                            f_del_m[pos] += 1
                            #still mismatch for 
                            f_ins_mm[pos] += 1
                    current_pos_ref += length #manuel
                else:
                    print(code, length, file=sys.stderr)
                    
    return (fm, rm, fmm, rmm, 
            f_ins_m, r_ins_m, f_ins_mm, r_ins_mm, 
            f_del_m, r_del_m, f_del_mm, r_del_mm)


def add_n(qgram_annotate, n, q):
    """Extend qgram_annotate by adding q-grams which contain Ns."""
    to_add = {}
    
    if n == 0:
        print('No q-grams with Ns to add' , file = sys.stderr)
        return qgram_annotate

    i = 0 #counter for status info    
    #consider each possible q-gram for the given length q and number n
    for qgram_with_n in _get_all_qgrams(['A','C','G','T','N'], ['A','C','G','T','N'], q, 1, n):
        qgram_with_n = qgram_with_n[0] #q-gram that contain at least one and at most n Ns
        i += 1
        if i % 20000000 == 0: 
            print(' %s q-grams with N considered' %(i), file = sys.stderr)

        #compute all concrete q-grams of n_qgram
        possible_qgrams = get_qgramlist(qgram_with_n)
        sb_table = [0,0,0,0] #initialize strand bias table
        for p_qgram in possible_qgrams:
            if qgram_annotate.has_key(p_qgram):
                sb_table = _add_listelements(sb_table, qgram_annotate[p_qgram][:])
        
        #does n containing q-gram corresponds to a combosed strand bias table?
        if sb_table != [0,0,0,0]: 
            to_add[qgram_with_n] = sb_table

    #extend qgram_annotate
    qgram_annotate.update(to_add)
    
    return qgram_annotate


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
            print(" %s / %s Strand Bias Scores calculated" %(i, len(qgram_annotate.keys())))
        
        #get p-value for the strand bias table of q-gram k
        p_value = get_pvalue(qgram_annotate[k][0], qgram_annotate[k][1], qgram_annotate[k][2], qgram_annotate[k][3])
        
        #compute negative logarithm (base 10) of p-value or, if necessary, set to maxint
        strand_bias_score = sys.maxint if p_value < 1/10.0**300 else -math.log(p_value, 10)
        
        results.append((k, qgram_annotate[k][0], qgram_annotate[k][1], qgram_annotate[k][2], qgram_annotate[k][3], strand_bias_score))
    
    return results


def output(results, genome, task):
    """Output the results"""
    print("Time before output: ",str(time.time()))
    print("###############################",task,"#########################")
    print("#Sequence", "Occurrence", "Forward Match", "Backward Match", "Forward Mismatch", "Backward Mismatch", "Strand Bias Score", "FER (Forward Error Rate)", 
            "RER (Reverse Error Rate), ERD (Error rate Difference)", sep = '\t')
    
    for seq, f_m, r_m, f_mm, r_mm, sb_score, fer, rer, erd in results:
        occ = count(seq, genome)
        print(seq, occ, f_m, r_m, f_mm, r_mm, sb_score, fer, rer, erd, sep = '\t')
    print("Time after output: ",str(time.time()))


def get_motifspace_size(q,n):
    """return length of search space according to equation which is mentioned in Section 3.1 of the paper"""
    return reduce(lambda x, y: x + (int(sc.comb(q, y, exact=True)) * 4**(q-y)), 
            [i for i in range(1, n+1)], int(sc.comb(q,0,exact=True)) * 4**(q-0))

def count(qgram, genome):
    """Count number of q-grams and its reverse complement in genome"""
    rev = reverse_complement(qgram)
    rev = rev.replace('N', '.')
    qgram = qgram.replace('N', '.')
    
    return  len([m.start() for m in re.finditer(r'(?=(%s))' %qgram, genome)] 
            + [m.start() for m in re.finditer(r'(?=(%s))' %rev, genome)])

def log(results, s, genome):
    """Result is a list [(seq, fm, rm, fmm, rmm, p_value)] """
    if not DEBUG_INDEL:
        return
    if s == "deletions":
        with open(LOG_DEL, 'w') as f:
            print("#Sequence", "Occurrence", "Forward Match", "Backward Match", 
                  "Forward Mismatch", "Backward Mismatch", "Strand Bias Score",
                  "FER","RER", sep = '\t', file = f)
            for r in results:
                seq = r[0]
                fer = float(r[3]) / (r[1] + r[3]) #forward error rate
                rer = float(r[4]) / (r[2] + r[4]) #forward error rate
                print(seq,count(seq,genome),str(r[1]),str(r[2]),str(r[3]),
                        str(r[4]),str(r[5]),str(fer),str(rer),sep='\t',file=f)


def ident(genome, genome_annotate, q, n, alpha=0.05, epsilon=0.03, delta=0.05):
    """Identify critical <q>-grams (with <n> Ns) with reference to significance and error rate""" 
    motifspacesize_log = math.log(get_motifspace_size(q, n), 10)
    alpha_log = math.log(float(alpha), 10)
    
    qgram_annotate, qgram_ins_annotate, qgram_del_annotate = get_annotate_qgram(
            genome, genome_annotate, q) #annotate each q-gram with Strand Bias Table
    add_n(qgram_annotate, n, q) #extend set of q-grams with q-grams containing Ns
    print("Number of in-, del keys before all_results [should match!]: ",
            len(qgram_annotate), len(qgram_ins_annotate), len(qgram_del_annotate))
    #annotate each q-gram with Strand Bias Score
    all_results = (get_sb_score(qgram_annotate),
            get_sb_score(qgram_ins_annotate),
            get_sb_score(qgram_del_annotate))
    print("Number of in-, del keys after all_results [should match!]: ",
            len(all_results[0]),len(all_results[1]), len(all_results[2]))
    log(all_results[2], "deletions", genome)
    #filter statistically significant motifs (Bonferroni Correction)
    all_sig_results = (filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results[0]),
            filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results[1]), 
            filter(lambda x: x[5] > motifspacesize_log - alpha_log, all_results[2])) 

    tasks = { 0 : 'snps', 1 : 'insertion', 2 : 'deletions' }
    for index, sig_results in enumerate(all_sig_results):
        #filter motifs with regards to background error rate (epsilon) error rate difference (delta)
        results = []
        for seq, f_m, r_m, f_mm, r_mm, p_value_score in sig_results:
            try:
                fer = float(f_mm) / (f_m + f_mm) #forward error rate
            except ZeroDivisionError:
                print("ZeroDivisionError fer:",count(seq,genome),tasks[index],
                        seq,str(f_m),str(r_m),str(f_mm),str(r_mm),str(p_value_score))
                continue
            try:
                rer = float(r_mm) / (r_m + r_mm) #reverse error rate
            except ZeroDivisionError:
                print("ZeroDivisionError rer:",count(seq,genome),tasks[index],
                        seq,str(f_m),str(r_m),str(f_mm),str(r_mm),str(p_value_score))
                continue

            #TODO remove, debug output
            if index > 0 and DEBUG_INDEL is True:
                print("(",tasks[index],count(seq,genome),seq,str(f_m),str(r_m),
                        str(f_mm),str(r_mm),str(p_value_score),str(fer),str(rer),")")

            if rer < epsilon: #filter sequences with too high epsilon (background error rate)
                erd = fer - rer #error rate difference
                if erd >= delta: #filter sequences with too low delta (error rate difference cutoff)
                    results.append((seq, f_m, r_m, f_mm, r_mm, p_value_score, fer, rer, erd)) 
            
        results.sort(key=lambda x: x[8],reverse=True) #sort by erd (error rate difference)
        output(results, genome, tasks[index])


if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    
    parser.add_option("-a", dest="alpha", default=0.05, type="float", help="FWER (family-wise error rate) alpha, default: 0.05")
    parser.add_option("-e", dest="epsilon", default=0.03, type="float", help="background error rate cutoff epsilon, default: 0.03")
    parser.add_option("-d", dest="delta", default=0.05, type="float", help="error rate difference cutoff delta, default: 0.05")
    parser.add_option("-c", dest="learn_chrom", default="chr1", help="chromosome that is used to derive Context Specific Errors, default: chr1")
    parser.add_option("-v", dest="version", default=False, action="store_true", help="show script's version")
    parser.add_option("-s", dest="serialize", default="", help="load serialized object instead of creating contingency tables for genome positions")
    
    (options, args) = parser.parse_args()
    
    if options.version:
            version = "version \"0.32\""
            print("")
            print(version)
            sys.exit()
    
    if len(args) != 4:
        parser.error("Sorry, exactly four parameters are required.")  
    
    #map arguments
    refpath = args[0]
    bampath = args[1]
    q = int(args[2])
    n = int(args[3])

    genome, options.learn_chrom = get_genome(refpath, options.learn_chrom)
    if options.serialize == "":
        print("Annotating (parsing) genome...")
        genome_annotate = get_annotate_genome(genome, bampath, options.learn_chrom)
        #barni TODO remove
        #object serialization
        with open('../serialized/genome_annotate_bsubtilis_indels.pkl', 'wb') as outpkl:
            pickle.dump(genome_annotate, outpkl, pickle.HIGHEST_PROTOCOL)
    else:
        #load input from serialized file
        print("Loading serialized object instead of annotating (parsing) genome...")
        with open(options.serialize, 'rb') as inpkl:
            genome_annotate = pickle.load(inpkl)
   # print(genome_annotate[4][165747],
   #        genome_annotate[5][165747],
   #        genome_annotate[6][165747],
   #        genome_annotate[7][165747])
   # print(genome_annotate[8][165747],
   #        genome_annotate[9][165747],
   #        genome_annotate[10][165747],
   #        genome_annotate[11][165747])


    ident(genome, genome_annotate, q, n, options.alpha, options.epsilon, options.delta)
