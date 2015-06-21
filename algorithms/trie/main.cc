#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include "dictionary.hh"
using namespace std;
int N = -1;

/**
 * @brief Extracts the first whole word from a line, starting at position pos
 * and using spaces/tabs as delimiters
 * @param line a set of words, space or tab delimited
 * @param pos the position in the line where to start looking for words
 * @return the found word, or an empty string
 */
string getNextWord(const string &line, unsigned &pos) {
    string word = "";
    unsigned j;
    for (j = pos; (j < line.length()) && (line[j] != ' ') && (line[j] != '\t'); j++) {
        word += line[j];
    }
    pos = j+1;
    return word;
}

/**
 * @brief Extracts the first whole word from a line, using spaces/tabs as delimiters
 * @param line a set of words, space or tab delimited
 * @return the found word, or an empty string
 */
string getNextWord(const string &line) {
    unsigned pos = 0;
    return getNextWord(line, pos);
}

void readVocabulary(const string &file, Dict &voc) {
    ifstream myfile (file);
    string line, word;
    int cnt = 0;

    if (myfile.is_open()) {
	    cout << "Loading file: " << file << endl;

        while(getline (myfile,line) && (cnt++ < N || N == -1))  {
            // skip comments and lines starting with space,
            // as they will be in most data files
            if (line[0] == '#' || line[0] == ' ')
                continue;

            // must extract the first whole word from the line
            word = getNextWord(line);

            voc.insert(word);
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
