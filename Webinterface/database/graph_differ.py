##
#
# Provide methods to plot different motifs between two platform technologies
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

def draw_differ(tech1, tech2, sequence):
    
    gene = tech1.gene  
    seg_loc = gt.seq_pos(gene, sequence)  
    
    tech1.sequence = sequence
    tech2.sequence = sequence
    
    'Init variables, arrays, trie'
    global gTrie1
    global gTrie2
    
    gTrie1 = trie()
    gTrie2 = trie()
    
    "Read a .csv file with motifs, build a Trie from it and store in gTrie"
    length = 0
    with open(tech1.motif_file.path) as f:
        motifs = f.readlines()
        for m in motifs[1:]:
            m = m.split(',')[0]
            length = len(m)
            gTrie1[m] = 1
    with open(tech2.motif_file.path) as f:
        motifs = f.readlines()
        for m in motifs[1:]:
            m = m.split(',')[0]
            length = len(m)
            gTrie2[m] = 1
            
    cMotifs = gt.getMatchingMotifs(gTrie1, tech2.motif_file.path)
    
    motifs = getDifferencMotifs(gTrie1, gTrie2, tech1, tech2, cMotifs)
    
    dMotifs = diff_motif_seq(sequence, motifs, length)
    dPos = diff_pos(sequence, motifs, length)
    
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
    yPos += 1
    'Draw Motifs'
    for i, m in enumerate(dMotifs):
        techNr = checkFindIn(gTrie1, gTrie2, m)
        if gt.checkShift(dMotifs[i], dPos[i], dPos[i-1]):
            yPos += 2 
        else:
            yPos = 1
        if techNr == 1:
            draw_motif(seg_loc,sequence, dMotifs[i], dPos[i], yPos, 'black')
            erd = gt.get_MotifERD(tech1.motif_file.path, dMotifs[i])
            print erd
            gt.draw_Information(ax, erd, dPos[i], yPos)
        if techNr == 2:
            draw_motif(seg_loc,sequence, dMotifs[i], dPos[i], yPos, 'grey')
            erd = gt.get_MotifERD(tech2.motif_file.path, dMotifs[i])
            print erd
            gt.draw_Information(ax, erd, dPos[i], yPos)

    'Read SNP'
    snps = gt.readSNP(tech1.snp_file.path)
    
    # "Check for each SNP if it belongs to a sequence and eventually plot"
    gt.process_SNP(tech1, sequence, seg_loc, dMotifs, dPos, snps, ax)
    
    'Last figure modifications'
    plt.ylabel("Difference")
    plt.legend(loc='upper right')
    plt.xticks(())
    plt.yticks(())
    
    save_folder = 'media/graphics/Bild_differ' + str(tech1.pk) + str(tech2.pk) + '.png'
    plt.savefig(save_folder)
    
    return save_folder
##
# Get a set M_a of all motifs
# Get a set M_c of all common motifs between different technologies
# M_d = M_a - M_c where M_d is set of all different motifs between technologies
##
def getDifferencMotifs(gTrie1, gTrie2, tech1, tech2, cMotifs):
          
    res = gt.getDifferentMotifs(gTrie1, tech2.motif_file.path)
    res.append(gt.getDifferentMotifs(gTrie2, tech1.motif_file.path))
    
    return res


def diff_pos(sequence, motifs, length):
    res = [ ]
        # iterate over ref string
    for i in range(0, len(sequence)-length+1):
        m = sequence[i:i+length]
        # if key is already there, append position
        if m in motifs:
            res.append(i)

    return res

'Retrieve common motifs which are in sequence'
def diff_motif_seq(sequence, motifs, length):
    res = [ ]
        # iterate over ref string
    for i in range(0, len(sequence)-length+1):
        m = sequence[i:i+length]
        # if key is already there, append position
        if m in motifs:
            res.append(m)

    return res

def checkFindIn(gTrie1, gTrie2, motif):
    
    if gTrie1.has_key(motif):
        return 1
    if gTrie2.has_key(motif):
        return 2
    
'Draw motif'
def draw_motif(sequence_locations, sequence, motif, xPos, yPos, color):
    buffer = 0.2
    
    for s in sequence_locations:
        plt.barh(yPos, width = len(motif), height=1, left=buffer+xPos, color=color)
        
        buffer += 10 + len(sequence)