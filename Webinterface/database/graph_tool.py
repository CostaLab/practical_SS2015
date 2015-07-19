##
# Provides basic methods to plot graphs
##
import csv

'Biopython libraries'
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq

'Graphical Libraries'
import matplotlib.pyplot as plt
import numpy as np
from units import unit

'Trie libraries and packages'
from optparse import OptionParser
from pytrie import SortedStringTrie as trie

'Import external python scripts'
import vcf

'Dictionary to determine color of Nucleobase'
NUC_TO_COLOR = {
    "A": "red",
    "a": "red",
    "C": "blue",
    "c": "blue",
    "G": "green",
    "g": "green",
    "T": "orange",
    "t": "orange",
}

'Draw Seqquence'
def draw_sequence(sequence_locations, sequence,  ax, yPos):
    buffer = 0.2
    
    for s  in sequence_locations:
        pos = int(s) - 1
        label = 'Position: ' + str(pos)
        ax.text(buffer, -0.2, label)
        
        for j, c in enumerate(sequence):
            color = NUC_TO_COLOR[c]
            plt.barh(0, width = len(c), height = 1, left = buffer+j, color=color)
        
        buffer += 10 + len(sequence)

'Draw all motifs'
def draw_all_motifs(tech, seg_loc,sequence, motifs, motif_pos, yPos,ax):
    i = 0
    draw_motif(seg_loc, sequence, motifs[i], motif_pos[i], yPos)
    erd = get_MotifERD(tech.motif_file.path, motifs[i])
    draw_Information(ax, erd, motif_pos[i], yPos)
    i += 1
    while i < len(motifs):
        'Draw Motifs'
        if checkShift(motifs[i], motif_pos[i], motif_pos[i-1]):
            yPos += 2
            draw_motif(seg_loc, sequence, motifs[i], motif_pos[i], yPos)
            erd = get_MotifERD(tech.motif_file.path, motifs[i])
            draw_Information(ax, erd, motif_pos[i], yPos)
        else:
            yPos = 1
            draw_motif(seg_loc, sequence, motifs[i], motif_pos[i], yPos)
            erd = get_MotifERD(tech.motif_file.path, motifs[i])
            draw_Information(ax, erd, motif_pos[i], yPos)
        i += 1
           
'Draw motif'
def draw_motif(sequence_locations, sequence, motif, xPos, yPos):
    buffer = 0.2
    
    for s in sequence_locations:
        plt.barh(yPos, width = len(motif), height=1, left=buffer+xPos, color='grey')
        
        buffer += 10 + len(sequence)
        
'Highlight SNP'
def draw_SNP(sequence_locations, sequence, snpALT, xPos, yPos):
    xBuffer = 0.2
    for s in sequence_locations:
        color = NUC_TO_COLOR[snpALT]
        plt.barh(yPos, width = len(snpALT), height=1, left=xBuffer+xPos, color=color)
        xBuffer += 10 + len(sequence)
        
'Draw annotation'
def draw_annotate(ax, information, xcoord, ycoord):
    ax.annotate(information, xy=(xcoord+0.5, ycoord + 1), xytext=(xcoord + 1, ycoord + 1),
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"),
                )
    
def draw_Information(ax, information, xcoord, ycoord):
    output = float("{0:.5f}".format(float(information)))
    ax.text(xcoord, ycoord+1.5, output,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=45
    )  
     
'Filter Error-Rate Difference of Motif'
def get_MotifERD(csv_file, motif):
    with open(csv_file) as f:
        for row in csv.reader(f, delimiter=','):
            if row[0] == motif:
                return row[3]
            
def trieGetMotifs(gTrie, s, q):
    """s is the reference (sub)string, q is length of motifs.
        Builds a trie where values are lists with positions of motif.
        """
#    global gTrie
    res = [ ]
    # iterate over ref string
    for i in range(0, len(s)-q+1):
        m = s[i:i+q]
        # if key is already there, append position
        if gTrie.has_key(m):
            res.append(m)

    return res

"Returns a list with motifs from file[fname_test] found in gTrie"
def getMatchingMotifs(gTrie, fname_test):
    
    res = []

    with open(fname_test) as f:
        motifs = f.readlines()        
        for m in motifs[1:]:
            m = m.split(',')[0]
            if gTrie.has_key(m):
                res.append(m)
            
                
    return res

"Returns a list with motifs from file[fname_test] found in gTrie"
def getDifferentMotifs(gTrie, fname_test):
    
    res = []

    with open(fname_test) as f:
        motifs = f.readlines()        
        for m in motifs[1:]:
            m = m.split(',')[0]
            if not gTrie.has_key(m):
                res.append(m)
                           
    return res

'Get relative position of a motif to a given sequence'
def getMotifPosition(gTrie, sequence, q):
    res = [ ]
    # iterate over ref string
    for i in range(0, len(sequence)-q+1):
        m = sequence[i:i+q]
        # if key is already there, append position
        if gTrie.has_key(m):
            res.append(i)

    return res

'Ermittle alle moeglichen Position der Sequence'
def seq_pos(gene, sequence):
    refFile = gene.refFile.path
    res = [ ]
    
    seq_record = SeqIO.read(refFile, "fasta")
    str_reference = seq_record.seq.rstrip()
    
    for i in range(0, len(str_reference)-len(sequence)+1):
        if str_reference[i:i+len(sequence)] == sequence:
            res.append(str(i+1))

    return res

'Define empty array to store possible snp variants'
def readSNP(vcf_file):
    res = [ ]
    res = vcf.VCFReader(open(vcf_file, 'rb'))
    
    return res

def process_SNP(tech, sequence, seg_loc, motifs, motif_pos, snps, ax):   

    for record in snps:
        annotate = False
        information = ""
        rel_POS = 0
        ySNP = 0
        
        for sl in seg_loc:
            min = int(sl)
            max = int(sl) + len(sequence)
            
            "Check position"
            if record.POS >= min and record.POS < max:
                "Check for motif"
                rel_POS = record.POS - min
                
                for i, m in enumerate(motif_pos):
                    min = int(m)
                    max = min + len(motifs[i])
                    if min <= rel_POS and rel_POS < max:
                        if annotate == False:
                            information += "REF: " + record.REF + "\n"
                            information += "ALT: " + prettyPrintALT(record.ALT) + "\n"
                            annotate = True
                        ySNP = checkSNP_pos(motifs[i], motif_pos, i)
                        draw_SNP(seg_loc, sequence, prettyPrintALT(record.ALT), rel_POS, ySNP)
 #                       erd = get_MotifERD(tech.motif_file.path, motifs[i])
 #                       if annotate == True:
 #                           information += motifs[i] + ": " + str(erd) + "\n"
        if annotate:
            draw_annotate(ax, information, rel_POS, ySNP)
            
# Prettty-Printer:  Input-> ['G'] Output G
def prettyPrintALT(alt):
    
    for a in alt:
        return a
    
    return alt

def buildTrie(tech, gTrie):
    "Read a .csv file with motifs, build a Trie from it and store in gTrie"
    with open(tech.motif_file.path) as f:
        motifs = f.readlines()
        for m in motifs[1:]:
            m = m.split(',')[0]
            length = len(m)
            gTrie[m] = 1
            
 
'Return true if previous motif is not overlapping with current motif'           
def checkShift(motif, m_pos, prev_pos): 
    length = len(motif)
    
    if m_pos == 0:
        return False
    if m_pos <= prev_pos + length:
        return True
    else:
        return False
    
'yPos of SNP has to be checked and determined'    
def checkSNP_pos(motif, m_pos, i):
    length = len(motif)
    yPos = 1
    
    while i > 0:
        if checkShift(motif, m_pos[i], m_pos[i-1]):
            yPos += 2
        else:
            return yPos
        i -= 1
        
    return yPos