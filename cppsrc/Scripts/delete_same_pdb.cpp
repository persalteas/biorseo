#include <iostream>
#include <sstream>
#include <fstream>
#include "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/cppsrc/json.hpp"
#include <typeinfo>
#include <set>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>

using namespace std;
using json = nlohmann::json;

void delete_redundant_pdb(const string& jsonlibrary, const string& fasta, const string& jsonoutfile) {
    std::ifstream lib(jsonlibrary);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;
    json js = json::parse(lib);

    std::ifstream file(fasta);
    string pdb, seq;
    std::getline(file, pdb);
    std::getline(file, seq);
    
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        vector<string> list_pdbs;
        vector<string> list_pdbs2;
        bool is_added = true;

        //cout << "id: " << id << endl;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            string test = it2.key();   

            if (!test.compare("pdb")) {
                vector<string> tab = it2.value();
                list_pdbs = tab;
                /*set<set<string>>::iterator iit;
                set<string>::iterator iit2;
                for(iit = list_pfams.begin(); iit != list_pfams.end(); iit++) {
                    for (iit2 = iit->begin(); iit2 != iit->end(); ++iit2) {
                        cout << *iit2 << endl;
                    }
                    cout << endl << endl;
                }*/
            } else {
                new_id[test] = it2.value();
            }  
        }
        
        if (count(list_pdbs.begin(), list_pdbs.end(), pdb.substr(6,pdb.size()))) {
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
    string jsonlibrary = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_final.json";
    string fasta = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/fasta/";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final.json";
    fasta = fasta + argv[1];
    delete_redundant_pdb(jsonlibrary, fasta, out);
    return 0;
}
    
