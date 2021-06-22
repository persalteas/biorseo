#include <iostream>
#include <sstream>
#include <fstream>
#include "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/cppsrc/json.hpp"
#include <typeinfo>
#include <set>
#include <algorithm>
#include <cstdio>
#include <vector>

using namespace std;
using json = nlohmann::json;

size_t count_delimiter(string& seq) {
    size_t count = 0;
    for(uint i = 0; i < seq.size(); i++) {
        char c = seq.at(i);
        if (c == '&') {
            count++;
        }
    }
    return count;
}

void add_delimiter(const string& jsonfile, const string& jsonoutfile) {
    std::ifstream lib(jsonfile);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;

    json js = json::parse(lib);
    
    //the list of pfam lists of the motif we want to count the inclusion in other motif
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        string test;
        string sequence;
        string contacts;
        bool is_change = false;

        //cout << "id: " << id << endl;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            test = it2.key();

            if (!test.compare("sequence")) {
                //cout << "sequence: " << it2.value() << endl;
                sequence = it2.value();
                new_id[test] = it2.value();
  
            } else if (!test.compare("contacts") ) {
                contacts = it2.value();
            } else {
                new_id[test] = it2.value();
            }  
        }
        string tmp = "";
        if (count_delimiter(contacts) != count_delimiter(sequence) && contacts.size() == sequence.size()) {
            for (uint i = 0; i < sequence.size(); i++) {
                if (sequence.at(i) == '&') {
                    tmp += "&";
                } else {
                    tmp += contacts.at(i);
                }
            }
        } else {
            tmp = contacts;
        }
        new_id["contacts"] = tmp;
        new_motif[id] = new_id;
        new_id.clear();
    }
    outfile << new_motif.dump(4) << endl;
    outfile.close();
    
}

int main()
{
    //183
    //cout << "------------------BEGIN-----------------" << endl;
    string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_06-06-2021.json";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_tmp.json";
    add_delimiter(jsonfile, out);

    //cout << "------------------END-----------------" << endl;
    return 0;
}
    
