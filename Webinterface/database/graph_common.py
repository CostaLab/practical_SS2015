##
#
# Provide methods to plot common motifs between two different technologies
#
#
##
import sys
import re
import csv

'Graphical Libraries'
import matplotlib.pyplot as plt
import numpy as np
from units import unit

'Trie libraries and packages'
from optparse import OptionParser
from pytrie import SortedStringTrie as trie

'Import external python scripts'
import graph_tool as gt

def draw_common(tech1, tech2, sequence):
    
    gene = tech1.gene  
    seg_loc = gt.seq_pos(gene, sequence)  
    
    tech1.sequence = sequence
    tech2.sequence = sequence
    
    'Init variables, arrays, trie'
    global gTrie
    gTrie = trie()
    
    'Prepare plot'
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5, forward=True)
    ax.invert_yaxis()

    'Prepare axis'
    ax.xaxis.tick_top()
    ax.xaxis.set_ticks_position('both') # THIS IS THE ONLY CHANGE
    ax.autoscale_view()

    yPos = 0
    'Draw the sequence first'
    gt.draw_sequence(seg_loc,sequence, ax, yPos)
    
    length = 0
    "Read a .csv file with motifs, build a Trie from it and store in gTrie"
    with open(tech1.motif_file.path) as f:
        motifs = f.readlines()
        for m in motifs[1:]:
            m = m.split(',')[0]
            length = len(m)
            gTrie[m] = 1
        
    motifs = gt.getMatchingMotifs(gTrie, tech2.motif_file.path)
    
    common_m = common_motif_seq(sequence, motifs, length)
    common_p = common_pos(sequence, common_m, length)
    
    yPos += 1
    'Draw Motifs'
    gt.draw_all_motifs(tech1, seg_loc,sequence, common_m, common_p, yPos, ax)
        
    'Read SNP'
    if not tech1.snp_file.path == None:
        snps1 = gt.readSNP(tech1.snp_file.path)
        gt.process_SNP(tech1, sequence, seg_loc, common_m, common_p, snps1, ax)
    if not tech2.snp_file.path == None:   
        snps2 = gt.readSNP(tech2.snp_file.path)
        gt.process_SNP(tech2, sequence, seg_loc, common_m, common_p, snps2, ax)
    
    
    
    'Last figure modifications'
    plt.ylabel("Common Motifs")
    plt.legend(loc='upper right')
    plt.xticks(())
    plt.yticks(())
    
    save_folder = 'media/graphics/Bild_common' + str(tech1.pk) + str(tech2.pk) + '.png'
    plt.savefig(save_folder)
    
    return save_folder

def common_pos(sequence, motifs, length):
    res = [ ]
        # iterate over ref string
    for i in range(0, len(sequence)-length+1):
        m = sequence[i:i+length]
        # if key is already there, append position
        if m in motifs:
            res.append(i)

    return res

'Retrieve common motifs which are in sequence'
def common_motif_seq(sequence, motifs, length):
    res = [ ]
        # iterate over ref string
    for i in range(0, len(sequence)-length+1):
        m = sequence[i:i+length]
        # if key is already there, append position
        if m in motifs:
            res.append(m)

    return res
    