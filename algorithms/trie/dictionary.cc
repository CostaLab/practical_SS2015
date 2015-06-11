#include "dictionary.hh"

Dict::~Dict() {
    for (auto it : this->getAllWords()) {
        trie.del(&trie, &trie, it.c_str());
    }
}

void Dict::insert(std::string w) {
//    std::cout << "inserting " << w << std::endl;
    Trie *t = trie.add(&trie, indexMap.size() + 1, w.c_str());
    if (t) {
        indexMap.insert(std::pair<unsigned, Trie *>(cnt++, t));
    }
}

void Dict::insert(unsigned nr, std::string w) {
    Trie *t = trie.add(&trie, nr, w.c_str());
    if (t) {
        indexMap.insert(std::pair<unsigned, Trie *>(cnt++, t));
    }
}

std::string Dict::find(unsigned index) {
    auto it = indexMap.find(index - 1);
    if (it != indexMap.end()) {
        return trie.find(it->second);
    }
    return std::string();
}

bool Dict::find(const std::string &w) {
//    getFrequency(w);
    return true;
}

// TODO implement this by creating a nice iterator
std::vector<std::string> Dict::getAllWords() {
    std::string s;
    std::vector<std::string> vec;
    for (unsigned i = 1; ; i++) {
        s = this->find(i);
        if (s.empty()) {
            break;   
        }
        vec.push_back(s);
    }
    return vec;
}

std::vector<std::string> Dict::getFrequency(std::string w) {
    std::vector<std::string> vec;
    trie.getFrequency(&trie, w.c_str(), "", vec);
    return vec;
}
