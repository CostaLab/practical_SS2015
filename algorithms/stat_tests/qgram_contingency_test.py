#!/usr/bin/python
import pysam as ps
import re
from optparse import OptionParser
from pprint import pprint

#file_fasta  = "/mnt/code/rwth/lab-bioinf/costalab/algorithms/data/bordetella_pertussis/hseq/NC_018518.1.fasta"
#file_bam  = "/mnt/code/rwth/lab-bioinf/costalab/algorithms/data/bordetella_pertussis/hseq/ERR142615.bam"

def reverse_complement(s, rev=True):
    """Return the reverse complement of a DNA sequence s"""
    _complement = dict(A="T", T="A", C="G", G="C", N="N")
    t = reversed(s) if rev else s
    try:
        rc = (_complement[x] for x in t)  # lazy generator expression
    except KeyError:
        return s
   
    return "".join(rc)  # join the elements into a string

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

def calc_contingency_table(file_fasta, file_bam, qgram_with_n):
    global ps
    fasta_file  = ps.Fastafile(file_fasta) # Opening fasta file
    sam_file    = ps.AlignmentFile(file_bam,"rb") # Opening AlignmentFile
    ref         = fasta_file.references[0]
    #qgram_with_n= "NCGTNAGAAC"

    qgrams      = _get_qgramlist_help(qgram_with_n, []) 
    qgrams_compl= [reverse_complement(x) for x in qgrams]

#start and end are according to the appeareance of motifs 
#(they don't appear before or after)
    start       = 0#2131000
    end         = 4023999#fasta_file.get_reference_length(ref) 
    total_end   = fasta_file.get_reference_length(ref) 
    ref_bases   = (fasta_file.fetch(ref, 0, total_end)).upper()
    poslist     = []
    poslist_rev = []
    snppos      = [] # position to be checked
    snppos_rev  = [] # position to be checked
#dictionary where position to genome is key and qgram is value
    pos_to_qgram_rev  = {}
    qgram_end_base = 'C'

    forward_match   = 0
    reverse_match   = 0
    forward_mismatch= 0
    reverse_mismatch= 0

##############################################################################
#populate lists forward
    for q in qgrams:
        l = [m.start() for m in re.finditer(r'(?=(%s))' % q, ref_bases)]
        for x in l:
            poslist.append(x) 
            snppos.append(x + 9) # check base at the end of motif (9 cause length is 10)
#populate lists for complemented
    for q in qgrams_compl:
        l = [m.start() for m in re.finditer(r'(?=(%s))' % q, ref_bases)]
        for x in l:
            poslist_rev.append(x) 
            snppos_rev.append(x) 
            pos_to_qgram_rev[x] = q

    print "--- qgrams to consider (forward)"
    pprint(qgrams)
    print "Positions of NCGTNAGAAC qgram:", poslist
    print "Check bases (pileups) at positions for NCGTNAGAAC qgram:", snppos

    print "--- qgrams (reverse complemented) to consider"
    print "Positions of NCGTNAGAAC qgram (complemented):", poslist_rev
    print "Check bases (pileups) at positions for NCGTNAGAAC qgram (complemented):", snppos_rev
    pprint(qgrams_compl)

##############################################################################
# Iterating through pileup
    print "Iterating through pileups..."
#for pileupcolumn in sam_file.pileup(ref, start, end): 
    for pileupcolumn in sam_file.pileup(None, None, None): 

        # forward (not complemented stuff)
        if (pileupcolumn.pos) in snppos:
            print "Pileup at position: %s, length: %s" % (pileupcolumn.pos, len(pileupcolumn.pileups))
            for ps in pileupcolumn.pileups:
                if not ps.is_del and not ps.is_refskip:
                    read_base = (ps.alignment.query_sequence[ps.query_position]).upper()
                else:
                    continue

                #print len(pileupcolumn.pileups), cnt
                if ps.alignment.is_reverse:
                    if read_base == qgram_end_base:
                        reverse_match += 1
                    else:
                        reverse_mismatch += 1
                else:
                    if read_base == qgram_end_base:
                        forward_match += 1
                    else:
                        forward_mismatch += 1

        # reverse complemented stuff
        elif (pileupcolumn.pos) in snppos_rev:
            print "Pileup at position (complemented): %s, length: %s" % (pileupcolumn.pos, len(pileupcolumn.pileups))
            for ps in pileupcolumn.pileups:
                if not ps.is_del and not ps.is_refskip:
                    read_base = (ps.alignment.query_sequence[ps.query_position]).upper()
                else:
                    continue

                if ps.alignment.is_reverse:
                    if read_base == pos_to_qgram_rev[pileupcolumn.pos][0]:
                        forward_match += 1
                    else:
                        forward_mismatch += 1
                else:
                    if read_base == pos_to_qgram_rev[pileupcolumn.pos][0]:
                        reverse_match += 1
                    else:
                        reverse_mismatch += 1


    print forward_match, reverse_match, forward_mismatch, reverse_mismatch
    fasta_file.close()
    sam_file.close()

if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    
    if len(args) != 3:
        parser.error("Sorry, exactly three parameter (fasta, bam, qgram with N) are required.")  

#    file_fasta      = args[0]
#    file_bam        = args[1]
#    qgram_with_n    = args[2]
    calc_contingency_table(args[0], args[1], args[2])
