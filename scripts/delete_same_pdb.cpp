#include <iostream>
#include <sstream>
#include <fstream>
#include "/local/local/BiorseoNath/cppsrc/json.hpp"
#include <typeinfo>
#include <set>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>

using namespace std;
using json = nlohmann::json;

/*
This script is use to create a new motif library without a motif that contains the same pdb as the sequence used in input for prediction
with BiORSEO.
*/
void delete_redundant_pdb(const string& jsonlibrary, const string& name, const string& jsonoutfile) {
    std::ifstream lib(jsonlibrary);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;
    json js = json::parse(lib);
    
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        vector<string> list_pdbs;
        bool is_added = true;

        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            string field = it2.key();   

            if (!field.compare("pdb")) {
                vector<string> tab = it2.value();
                list_pdbs = tab;
            } else {
                new_id[field] = it2.value();
            }  
        }

        if (count(list_pdbs.begin(), list_pdbs.end(), name.substr(0, name.size()-2))) {
            is_added = false;
        }
        if (is_added) {      
            new_id["pdb"] = list_pdbs;     
            new_motif[id] = new_id;
        }
        new_id.clear();
    }
    outfile << new_motif.dump(4) << endl;
    outfile.close(); 
}

int main(int argc, char** argv)
{
    string jsonlibrary = "/local/local/BiorseoNath/data/modules/ISAURE/motifs_final.json";
    string out = "/local/local/BiorseoNath/data/modules/ISAURE/bibliotheque_a_lire/motifs_final.json";
    string name = argv[1];
    delete_redundant_pdb(jsonlibrary, name, out);
    return 0;
}
    
