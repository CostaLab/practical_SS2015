#!/usr/bin/python
###############################################################################
# You can use this script in 2 way:                                           #
#   1) as a module: you can just use the two methods defined                  #  
#   2) as an exec script: run the script with 2 parameters                    #  
###############################################################################

from __future__ import print_function
from optparse import OptionParser
from pytrie import SortedStringTrie as trie
import sys

def getMatchingMotifs(fname_test):
    global gTrie

    with open(fname_test) as f:
        motifs = f.readlines()        
        for m in motifs[1:]:
            m = m.split(',')[0]
            if T.has_key(m):
                print(m)

def buildTrie(fname_ref):
    global gTrie

    with open(fname_ref) as f:
        motifs = f.readlines()        
        for m in motifs[1:]:
            m = m.split(',')[0]
            T[m] = 1
    

if __name__ == '__main__':
    parser          = OptionParser()
    (options, args) = parser.parse_args()
    
    if len(args) > 1:
        print("Two files provided, processing...", file=sys.stderr)

        fname_ref   = args[0]
        fname_test  = args[1]
        T           = trie()
        
        with open(fname_ref) as f:
            motifs = f.readlines()        
            for m in motifs[1:]:
                m = m.split(',')[0]
                T[m] = 1
                
        with open(fname_test) as f:
            motifs = f.readlines()        
            for m in motifs[1:]:
                m = m.split(',')[0]
                if T.has_key(m):
                    print(m)
    else:
        global gTrie 
        gTrie = trie()
