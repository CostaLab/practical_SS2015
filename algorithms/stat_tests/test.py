#!/usr/bin/python

a = ["AAA", "NAA", "NNT", "TAC", "TNC"]
b = ["AAA", "NNA", "CGC", "TAN", "TNN"]

#result should be AAA, NAA, NNA, NNT, TNN, TNC, TAC, TAN

#return 
def compare(s, w):
    for i in range(len(s)):
        if s[i] == w[i]:
            continue
        if s[i] == 'N':
            return -1
        else:
            return 1
    return 0

print sorted(a, cmp = compare)
print sorted(b, cmp = compare)
