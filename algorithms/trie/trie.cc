#include "trie.hh"

Trie::Trie():
    pred(0),
    first_indx(-1),
    occur(0),
    son(0),
    nr(0)
{
    bases = { 'A', 'C', 'T', 'G', 'N' };
    std::memset(next,0,sizeof(next));
    punctuation['?'] = OFFSET + 1;
    punctuation['!'] = OFFSET + 2;
    punctuation['['] = OFFSET + 3;
    punctuation[']'] = OFFSET + 4;
    punctuation['('] = OFFSET + 5;
    punctuation[')'] = OFFSET + 6;
    punctuation['\''] = OFFSET + 7;
    punctuation['-'] = OFFSET + 8; 
    punctuation[','] = OFFSET + 9; 
    punctuation['.'] = OFFSET + 10; 
    punctuation[';'] = OFFSET + 11; 
    punctuation[':'] = OFFSET + 12; 
    punctuation['*'] = OFFSET + 13; 
    punctuation['_'] = OFFSET + 14; 
    punctuation['`'] = OFFSET + 15; 
    punctuation['"'] = OFFSET + 16; 
} 

int Trie::getPos(const char *w) {
//    std::cout << "getpos " << w << std::endl;
    if (*w >= 'A' && *w <='Z') {
        return *w-'A'; 
    }
    else if (*w <= '9' && *w >= '0') {
        return OFFSET_NR + (*w - '0') + 1;
    }
    else return punctuation[*w];
}

int Trie::getPos(const char w) {
//    std::cout << "getpos ok " << w << std::endl;
    if (w >= 'A' && w <='Z') {
        return w-'A'; 
    }
    else if (w <= '9' && w >= '0') {
        return OFFSET_NR + (w - '0') + 1;
    }
    else return punctuation[w];
}
/*
    Returns a pointer to the leave if insertion was successful,
    and 0 otherwise (word already in the trie).
*/
Trie* Trie::add(Trie *t, int nr, const char *w) {
    if (*w=='\0') {
        if (t->occur) {
            t->occur += nr;
            t->nr ++;
            return 0;
        }
        t->occur += nr;
        t->nr ++;
        return t;
    }
    int pos;
    pos = getPos(w);
    if (!t->next[pos]) {
        ++t->son;
        t->next[pos] = new Trie;
        t->next[pos]->pred = t;
        t->next[pos]->c = *w;
    }
    return add(t->next[pos], nr, w+1);
}

std::string Trie::find(Trie *t) {
    if (!t->pred) {
        return std::string();
    }
    return (find(t->pred)) + t->c;
}

void Trie::getFrequency(Trie *t, const char *w, std::string const_w, std::vector<std::string> &result) {
    if (!t) {
        return; 
    }

//    std::cout << "w " << w << " const_w " << const_w << std::endl;

    if (*w=='\0') {
        if (t->occur > 0) {
            result.push_back(const_w);
        }
        return;
    }

    for (int i = 0; i < 5; ++i) {
        char c[2];
        sprintf(c,"%c",bases[i]);
        if (*w == bases[i] || *w == 'N' || bases[i] == 'N') {
            unsigned pos = getPos(bases[i]);
            getFrequency(t->next[pos], w+1, const_w + std::string((const char*) c), result);
        }
    }
}

int Trie::del(Trie *T, Trie *t, const char *w)
{
    int pos = getPos(w);
    if (*w == '\0')
        --t->nr;
    else if (del(T,t->next[pos],w+1)) {
        t->next[pos]=0;
        --t->son;
    }
    if (!t->son && !t->nr && t!=T) {
        delete t;
        return 1;
    }
    return 0;
}

