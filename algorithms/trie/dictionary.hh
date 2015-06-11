#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <iostream>
#include <unordered_map>
#include "trie.hh"
typedef std::unordered_map<unsigned, Trie *> IndexMap;

class Dict {
    public:
        Dict() : cnt(0) {}
        ~Dict();

        void        insert(std::string w);
        void        insert(unsigned nr, std::string w);
        std::string find(unsigned ind);
        bool        find(const std::string &w);
        std::vector<std::string>    getFrequency(std::string w);
        std::vector<std::string> getAllWords();
        
    private:
        Trie trie;
        IndexMap indexMap;
        unsigned long cnt;
};

#endif // DICTIONARY_H
