#!/usr/bin/python
import khmer_dict
import numpy
from timeit import timeit
import sys
import resource

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
                #not too many Ns
                if letter == 'N' and el.count('N') >= n:
                    continue
                #not too less Ns
                if length - level <= 1 and el.count('N') == 0 and letter != 'N':
                    continue

                for r in _get_all_qgrams(alphabet, [el + letter], length, level+1, n):
                    yield r
q = 10 
d = {} 
cnt = 0
for qgram_with_n in _get_all_qgrams(['A','C','G','T','N'], ['A','C','G','T','N'], q, 1, 2):
#for qgram_with_n in _get_all_qgrams(['A','C','G','T'], ['A','C','G','T'], q, 1, 1):
    d[qgram_with_n[0]] = [1,2,3,4]
    if cnt % 100000 == 0:
        print cnt, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
    cnt += 1

print memory_usage()
print "#qgrams that were processed:", cnt
