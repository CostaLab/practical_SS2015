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
int N = -1;

string getNextWord(const string &line, unsigned &pos) {
    string word = "";
    unsigned j;
    for (j = pos; (j < line.length()) && (line[j] != ' ') && (line[j] != '\t'); j++) {
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
	    cout << "Loading file: " << file << endl;

        while(getline (myfile,line) && (cnt++ < N || N == -1))  {
            pos = 0;
            while (pos < line.length()) {
                word = getNextWord(line, pos);
                voc.insert(word);
            }
        }

        cout << "\t-> Words loaded: " << cnt << endl;

        myfile.close();
    }
    else {
        cout << "Unable to open vocabulary file: " << file;
    }
}

int main (int argc, char *argv[]) {
    if (argc < 4)
    {
        cout << "Usage: " << argv[0] << " file1 file2 [fileN..] file_out\n";
        exit(-1);
    }

    std::vector<Dict> inputs;

    for (int i = 1; i < argc - 1; i++)
    {
        Dict temp_dict;
        readVocabulary(argv[i], temp_dict);
        inputs.push_back(temp_dict);
    }

    // load first dictionary in a set
    std::vector<std::string> first_dict_words = inputs[0].getAllWords();
    std::set<std::string> common(first_dict_words.begin(), first_dict_words.end());

    // we are going to need a temporary structure for
    // common words
    std::set<std::string> temp_common;

    // will keep the results from a single search in the trie
    std::vector<std::string> res;

    for (unsigned ii = 1; ii < inputs.size(); ii++)
    {
        for (auto word : common)
        {
            res = inputs[ii].getFrequency(word);

            if (res.size() > 0)
            {
                temp_common.insert(word);

                for (auto common_word : res)
                    temp_common.insert(common_word);
            }
        }

        common = temp_common;

        temp_common.clear();
    }

    // now we have the common words for all files,
    // let's write them to file

    cout << "Writing common words to file: " << argv[argc-1] << endl;

    ofstream outfile (argv[argc-1]);

    int i = 0;
    for (auto word : common)
    {
        outfile << word;
        outfile << std::endl;

        i++;
    }

    cout << "\t-> Words written: " << i << endl;

    return 0;
}
