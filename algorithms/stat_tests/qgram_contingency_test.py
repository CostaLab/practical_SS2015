#!/usr/bin/python
from __future__ import print_function
import pysam as ps
import re
import sys
from optparse import OptionParser
from pprint import pprint

#file_fasta  = "/mnt/code/rwth/lab-bioinf/costalab/algorithms/data/bordetella_pertussis/hseq/NC_018518.1.fasta"
#file_bam  = "/mnt/code/rwth/lab-bioinf/costalab/algorithms/data/bordetella_pertussis/hseq/ERR142615.bam"
verbose = False

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

    N           = len(qgram_with_n)
    qgrams      = _get_qgramlist_help(qgram_with_n, []) 
    qgrams_compl= [reverse_complement(x) for x in qgrams]

#start and end are according to the appeareance of motifs 
#(they don't appear before or after)
    start       = 0#2131000
    end         = 4023999#fasta_file.get_reference_length(ref) 
    total_end   = fasta_file.get_reference_length(ref) 
    genome      = (fasta_file.fetch(ref, 0, total_end)).upper()
    poslist     = []
    poslist_rev = []
    snppos      = [] # position to be checked (motif_start + motif_length - 1)
    snppos_rev  = [] # position to be checked
#dictionary where position to genome is key and qgram is value
    pos_to_qgram      = {}
    pos_to_qgram_rev  = {}

    forward_match   = 0
    reverse_match   = 0
    forward_mismatch= 0
    reverse_mismatch= 0

##############################################################################
#populate lists forward
    for q in qgrams:
        l = [m.start() for m in re.finditer(r'(?=(%s))' % q, genome)]
        for x in l:
            poslist.append(x) 
            snppos.append(x + N - 1) # check base at the end of motif (9 cause length is 10)
            pos_to_qgram[x] = q
#populate lists for complemented
    for q in qgrams_compl:
        l = [m.start() for m in re.finditer(r'(?=(%s))' % q, genome)]
        for x in l:
            poslist_rev.append(x) 
            snppos_rev.append(x) 
            pos_to_qgram_rev[x] = q

    if verbose:
        print("--- qgrams to consider (forward)")
        pprint(qgrams)
        print( "Positions of NCGTNAGAAC qgram:")
        pprint(poslist)
        print( "Check bases (pileups) at positions for NCGTNAGAAC qgram:")
        pprint(snppos)

        print( "--- qgrams (reverse complemented) to consider")
        print( "Positions of NCGTNAGAAC qgram (complemented):", str(poslist_rev))
        print( "Check bases (pileups) at positions for NCGTNAGAAC qgram (complemented):", str(snppos_rev))
        pprint(qgrams_compl)

##############################################################################
# Iterating through pileup
    print ("Iterating through pileups...", file=sys.stderr)
#for pileupcolumn in sam_file.pileup(ref, start, end): 
    for pileupcolumn in sam_file.pileup(None, None, None): 

        # forward (not complemented stuff)
        if (pileupcolumn.pos) in snppos: # this is motif_start + motif_length
            qgram_first_pos = pileupcolumn.pos - N + 1
            qgram_last_pos  = pileupcolumn.pos 
            q = pos_to_qgram[qgram_first_pos]
            if verbose:
                print("Pileup at position: %s, \n\tpileup length: %s \n\tmaps to qgram: %s \
                    \n\tqgram_start_pos: %s\n\tqgram_end_pos (==pileup_pos): %s" % (qgram_last_pos, 
                    len(pileupcolumn.pileups), q, qgram_first_pos, qgram_last_pos), file=sys.stderr)
            if (q != genome[qgram_first_pos : qgram_first_pos + len(q)]):
                print("Error. Qgram [%s] does not match position / substring in genome [%s]:" % 
                        (q, genome[pileupcolumn.pos : pileupcolumn.pos + len(q)]), file=sys.stderr)
                exit(-1)

            for ps in pileupcolumn.pileups:
                #TODO remove, debug
                if ps.alignment.cigar is None:
                    print("wtf cigar not available: %s", (str(ps.alignment)))
                # if read is unmapped, discard / continue
                if ps.alignment.is_unmapped:
                    continue  
                if not ps.is_del and not ps.is_refskip:
                    read_base = (ps.alignment.query_sequence[ps.query_position]).upper()
                else:
                    continue

                if ps.alignment.is_reverse:
                    if read_base == q[-1:]:
                        reverse_match += 1
                    else:
                        reverse_mismatch += 1
                else:
                    if read_base == q[-1:]:
                        forward_match += 1
                    else:
                        forward_mismatch += 1

        # reverse complemented stuff
        elif (pileupcolumn.pos) in snppos_rev:
            q = pos_to_qgram_rev[pileupcolumn.pos]
            if verbose:
                print("Pileup at position: %s (reverse)\n\tpileup length: %s \n\tmaps to qgram: %s\
                    \n\tqgram reverse complemented: %s" \
                    % (pileupcolumn.pos, len(pileupcolumn.pileups), reverse_complement(q), q),
                    file=sys.stderr)
            if (q != genome[pileupcolumn.pos : pileupcolumn.pos + len(q)]):
                print("Error reverse. Qgram [%s] does not match position / substring in genome [%s]:" % 
                        (q, genome[pileupcolumn.pos : pileupcolumn.pos + len(q)]), file=sys.stderr)
                exit(-1)

            for ps in pileupcolumn.pileups:
                #TODO remove, debug
                if ps.alignment.cigar is None:
                    print("wtf cigar not available: %s", (str(ps.alignment)))
                # if read is unmapped, discard / continue
                if ps.alignment.is_unmapped:
                    continue  

                if not ps.is_del and not ps.is_refskip:
                    read_base = (ps.alignment.query_sequence[ps.query_position]).upper()
                else:
                    continue

                if ps.alignment.is_reverse:
                    if read_base == q[0]:
                        forward_match += 1
                    else:
                        forward_mismatch += 1
                else:
                    if read_base == q[0]:
                        reverse_match += 1
                    else:
                        reverse_mismatch += 1


    print(str(forward_match), str(reverse_match), str(forward_mismatch), str(reverse_mismatch))
    fasta_file.close()
    sam_file.close()

if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()
    
    if len(args) != 3:
        parser.error("Sorry, exactly three parameter (fasta, bam, qgram with N) are required.")  
    print("Change verbosity in code if you want...")

    calc_contingency_table(args[0], args[1], args[2])
