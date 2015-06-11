#ifndef TRIE_H
#define TRIE_H

#include<map>
#include<vector>
#include<cstdio>
#include<cstring>
#include<string>
#include<iostream>
#define TRIE_SIZE 26 + 10 + 16
#define OFFSET 35
#define OFFSET_NR 25 

// Trie data structure
class Trie {
    public:
        Trie();
        ~Trie() {}

        Trie*   add(Trie *t, int ind, const char *c);
        std::string find(Trie *t);
        void getFrequency(Trie *t, const char *w,std::string const_w, std::vector<std::string> &result);
        int     del(Trie *T,Trie *t,const char *w);
        int     getPos(const char *w);
        int     getPos(const char w);
    private:

        Trie *next[TRIE_SIZE];
        Trie *pred;
        char c;
        int first_indx;
        int occur;
        int son;
        int nr;
        std::map<char, int> punctuation;
        std::vector<char> bases;
};
#endif // TRIE_H
