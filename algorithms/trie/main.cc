// FIXME We assume that the training / test data are sorted according 
// to class names

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include "dictionary.hh"
using namespace std;
int N;

string getNextWord(const string &line, unsigned &pos) {
    string word = "";
    unsigned j;
    for (j = pos; (j < line.length()) && (line[j] != ' '); j++) {
        word += line[j];
    }
    pos = j+1;
    return word;
}

void readVocabulary(const string &file, Dict &voc) {
    ifstream myfile (file);
    string line, word;
    int cnt = 0;
    unsigned pos;
    if (myfile.is_open()) {
        while(getline (myfile,line) && (cnt++ < N || N == -1))  {
            pos = 0;
            while (pos < line.length()) {
                word = getNextWord(line, pos);
                voc.insert(word);
            }
        }
        myfile.close();
    }
    else {
        cout << "Unable to open vocabulary file";
    }
}

int main (int argc, char *argv[]) {
    string line;
    string file1_name;
    string file2_name;
    string out_name;
//    string param;
    Dict dic1;
    Dict dic2;

    N = -1;
    if(argc > 1) {
//        param = argv[4];
        file1_name = argv[1];
        file2_name = argv[2];
        out_name = argv[3];
    }

    readVocabulary(file1_name,dic1);

//    cout << "length: " << res.size() << endl;
//    for (unsigned i = 0; i < res.size(); ++i) {
//        std::cout << res[i] << std::endl;
//    }
    ifstream myfile (file2_name);
    ofstream outfile (out_name);
    string word;
    int cnt = 0;
    unsigned pos;
    std::vector<std::string> res;
    if (myfile.is_open()) {
        while(getline (myfile,line) && (cnt++ < N || N == -1))  {
            pos = 0;
            while (pos < line.length()) {
                word    = getNextWord(line, pos);
                res     = dic1.getFrequency(word);
                if (res.size() > 0) {
                    outfile << (word);
                    outfile << (std::endl);
                    for (unsigned i = 0; i < res.size(); ++i) {
                        outfile << (res[i]);
                        outfile << ("\n");
                    }
                }
            }
        }
        myfile.close();
        outfile.close();
    }
    return 0;
}
