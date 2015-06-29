#!/usr/bin/python
import sys
import math
import resource
import gc
import jellyfish

def memory_usage():
    """Memory usage of the current process in kilobytes."""
    status = None
    result = {'peak': 0, 'rss': 0}
    try:
        status = open('/proc/self/status')
        for line in status:
            parts = line.split()
            key = parts[0][2:-1].lower()
            if key in result:
                result[key] = int(parts[1])
    finally:
        if status is not None:
            status.close()
    return result

def _get_all_qgrams(alphabet, erg, length, level, n):
    """Return all possible q-grams of the given length, over the given alphabet
    and with at least one and at most n Ns"""
    if length == level:
        yield erg
    else:
        for letter in alphabet:
            for el in erg:
                for r in _get_all_qgrams(alphabet, [el + letter], length, level+1, n):
                    yield r

q = 15
cnt = 0
k = jellyfish.HashCounter(int(math.pow(2,23)), q)
for qgram_with_n in _get_all_qgrams(['A','C','G','T'], ['A','C','G','T'], q, 1, 1):
    if not k.add(jellyfish.MerDNA(qgram_with_n[0]), 1):
        print "Error."
        exit()
    if cnt % 10000 == 0:
        print cnt, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
    cnt += 1
print memory_usage()
#k.update_add(jellyfish.MerDNA("AAATGCAATG"), 1)
#print k[jellyfish.MerDNA("AAATGCAATG")]
print "#qgrams in hash table:", cnt
