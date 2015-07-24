#!/usr/bin/python
import numpy 
import math
import sys
import scipy.misc as sc

class Khmer:
    def __init__(self, k, n):
        # l is log2(motifspace_size) + 1
        self.MAX_REPROBE = 62
        self.l      = int(math.log(self.get_motifspace_size(k,n),2)) + 1
        self.k      = k
        self.M      = pow(2,self.l)
        self.not_added = 0
        # HASH_TABLE
        # has 5 elements, table[0] being the key, rest is fm, rm, fmm, rmm
        self.table 	= numpy.empty(self.M, dtype=object)
        A           = numpy.random.random_integers(0,1,(3*k, 3*k))
        while numpy.linalg.det(A) == 0:
            A       = numpy.random.random_integers(0,1,(3*k, 3*k))
        Ainv   = numpy.linalg.inv(A)
        self.A_list = []
        # transform each row of A to an int, put them in a list
        for r in A:
            v = 0
            for bit in r:
                v = (v << 1) | bit
            self.A_list.append(v)

    def getBitArray(self, key):
        bitarray = 0
        for i, c in enumerate(key):
            bitarray = bitarray << 3 
            if c == 'A':
                bitarray += 0
            elif c == 'C':
                bitarray += 1
            elif c == 'G':
                bitarray += 2
            elif c == 'T':
                bitarray += 3
            elif c == 'N':
                bitarray += 4
        return bitarray

    
    def __len__(self):
        l = 0
        for v in self.table:
            if v != None:
                l += 1
        return l


    def __getitem__(self, string_key):
        # create list of bits (3 * k) from key
        bit_key          = self.getBitArray(string_key)
        # compute hash function as product of matrix and bitarray
        func = 0
        nr_columns = 3 * self.k - 1
        for ind, r in enumerate(self.A_list):
            bitand = bit_key & r
            func += ((bin(bitand).count('1')%2) << (nr_columns - ind))
            #print bit_key, r, bitand, func, nr_columns

        hash_higher_bits = func / self.M # take highest 3*k - l bits
        hash_func        = func % self.M # hash(m) = f(m) % M
        mod_key = (1<<(3*self.k-self.l))
        i = 0

        while True:
            if i >= self.MAX_REPROBE:
                return None
            hash_pos = (hash_func + i*(i+1)/2) % self.M # take lowest l bits of [func + reprobe(i)]
            # if pos is empty, key is not in table
            if self.table[hash_pos] is None:
                return None
            cur_key  = self.table[hash_pos][0] % (mod_key) 
            i       += 1
            if cur_key == hash_higher_bits:
                break

        return self.table[hash_pos][1:]


    def __setitem__(self, string_key, val):
        # create list of bits (3 * k) from key
        bit_key          = self.getBitArray(string_key)
        # compute hash function as product of matrix and bitarray
        func = 0
        nr_columns = 3 * self.k - 1
        for ind, r in enumerate(self.A_list):
            bitand = bit_key & r
            func += ((bin(bitand).count('1')%2) << (nr_columns - ind))
            #print bit_key, r, bitand, func, nr_columns

        hash_higher_bits = func / self.M # take highest 3*k - l bits
        hash_func        = func % self.M # hash(m) = f(m) % M

        i       = 0 # i is same as in the paper : reprobe(i)
        cur_key = 0
        mod_key = (1<<(3*self.k-self.l))

        while True:
            if i >= self.MAX_REPROBE:
                self.not_added += 1
                return False
            hash_pos = (hash_func + i*(i+1)/2) % self.M # take lowest l bits of [func + reprobe(i)]
            cur_key  = self.table[hash_pos][0] % (mod_key) if self.table[hash_pos] != None else None
            #print "cur_key", cur_key
            i       += 1
            if cur_key is None: #or cur_key == hash_higher_bits:
                break
        self.table[hash_pos] = [int(i << (3*self.k-self.l)) + hash_higher_bits] + val

#        print "f(m): ", bin(func)
#        print "hash(m) = f(m) % M: ", bin(hash_func)
#        print "Hash pos: ", bin(hash_pos)
#        print "pure bitkey: ", bit_key, bin(bit_key)
#        print "KEY: ", bin(self.table[hash_pos][0])
#        print "table: ", self.table[hash_pos], "\n\n"
        return True


    def get_motifspace_size(self, q,n):
        return reduce(lambda x, y: x + (int(sc.comb(q, y, exact=True)) * 4**(q-y)), 
                [i for i in range(1, n+1)], int(sc.comb(q,0,exact=True)) * 4**(q-0))
