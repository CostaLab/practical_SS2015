#!/bin/bash

import itertools

def load_file(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue

            # Motif, Occurrences, FM, BM, FMM, BMM, SBS, FER, RER, ERD
            fields = line.split()

            motif = fields[0]
            opts  = fields[1:6]

            # accumulate tables of same motif (there shouldn't be any here,
            # unless the file was a naive merge of multiple results files)
            if motif in d:
                d[motif].append(opts)
            else:
                d[motif] = [opts]

    for l in d:
        d[l] = [sum(sb) for sb in itertools.izip(*(d[l]))]

merge_dicts(list_of_dicts, d):
    if len(list_of_dicts) == 0:
        return d
    
    d2 = lists_of_dicts[0]
    
    d = {k:v.append(d2[k]) for k,v in d if k in d2}

    return merge_dicts(lists_of_dicts[1:], d)


