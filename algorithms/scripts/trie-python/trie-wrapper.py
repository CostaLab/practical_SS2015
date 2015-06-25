#!/usr/bin/python
###############################################################################
# You can use this script in 1 way:                                           #
#   1) as a module: you can just use the methods defined                      #  
#   //2) as an exec script: run the script with 2 parameters                  #  
###############################################################################

from __future__ import print_function
from optparse import OptionParser
from pytrie import SortedStringTrie as trie
import sys

def getMotifPositions(m):
    """Returns a list with positions of motif argument m"""
    global gTrie

    if gTrie.has_key(m):
        return gTrie[m]

    return []


def getMatchingMotifs(fname_test):
    """Returns a list with motifs from file[fname_test] found in gTrie"""
    global gTrie
    res = []

    with open(fname_test) as f:
        motifs = f.readlines()        
        for m in motifs[1:]:
            m = m.split(',')[0]
            if gTrie.has_key(m):
                res.append(m)
                #print(m)
    return res

def trieWithPosFromString(s, q):
    """s is the reference (sub)string, q is length of motifs.
    Builds a trie where values are lists with positions of motif.
    """
    global gTrie

    # iterate over ref string
    for i in range(len(s) - q + 1):
        m = s[i:q]
        # if key is already there, append position
        if gTrie.has_key(m):
            gTrie[m].append(i)
        else:
            # create a list containing only current position
            gTrie[m] = [i]

def trieFromFile(fname_ref):
    """Read a .csv file with motifs, build a Trie from it and store in gTrie"""
    global gTrie

    with open(fname_ref) as f:
        motifs = f.readlines()        
        # if first line already contains a motif, change 1 to 0 in line below
        for m in motifs[1:]:
            m = m.split(',')[0]
            gTrie[m] = 1
    

if __name__ == '__main__':
    """main() does nothing, except creating the global gTrie """
    parser          = OptionParser()
    (options, args) = parser.parse_args()

    global gTrie 
    gTrie = trie()
    
#    if len(args) > 1:
#        print("Two files provided, processing...", file=sys.stderr)
#
#        fname_ref   = args[0]
#        fname_test  = args[1]
#        T           = trie()
#        
#        with open(fname_ref) as f:
#            motifs = f.readlines()        
#            for m in motifs[1:]:
#                m = m.split(',')[0]
#                T[m] = 1
#                
#        with open(fname_test) as f:
#            motifs = f.readlines()        
#            for m in motifs[1:]:
#                m = m.split(',')[0]
#                if T.has_key(m):
#                    print(m)
#    else:
#        global gTrie 
#        gTrie = trie()
