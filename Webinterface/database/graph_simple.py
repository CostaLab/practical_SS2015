##
#
# Provide methods to plot motif_graph
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
import vcf
import graph_tool as gt

def draw_simple(tech, sequence):
        
    'Init variables, arrays, trie'
    global gTrie
    gTrie = trie()
    
    snps = [ ]

    gene = tech.gene  
    seg_loc = gt.seq_pos(gene, sequence)  
        
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
    with open(tech.motif_file.path) as f:
        motifs = f.readlines()
        for m in motifs[1:]:
            m = m.split(',')[0]
            length = len(m)
            gTrie[m] = 1

    'Get Motifs and their positions'
    motifs = gt.trieGetMotifs(gTrie, sequence, length)
    motif_pos = gt.getMotifPosition(gTrie, sequence, length)
    yPos += 1
    
    gt.draw_all_motifs(tech, seg_loc,sequence, motifs, motif_pos, yPos, ax)    
            
    'Read SNP'
    if not tech.snp_file == None:
        snps = gt.readSNP(tech.snp_file.path)
    
    "Check for each SNP if it belongs to a sequence and eventually plot"
    gt.process_SNP(tech, sequence, seg_loc, motifs, motif_pos, snps, ax)
    
    'Last figure modifucations'
    plt.ylabel("Motifs")
    plt.legend(loc='upper right')
    plt.xticks(())
    plt.yticks(())
    

    save_folder = 'media/graphics/Bild' + str(tech.pk) + '.png'
    plt.savefig(save_folder)
    
    return save_folder


