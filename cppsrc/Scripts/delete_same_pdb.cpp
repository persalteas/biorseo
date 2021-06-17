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

void delete_redundant_pdb(const string& jsonfile, const string& jsontest, const string& jsonoutfile) {
    std::ifstream lib(jsonfile);
    std::ifstream lib2(jsontest);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;
    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
    //the list of pfam lists of the motif we want to count the inclusion in other motif
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
        //cout << "-------begin---------" << endl;
        
        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();

            //cout << "id: " << id << " / id2: " << id2 << endl;
            for (auto it4 = js2[id2].begin(); it4 != js2[id2].end(); ++it4) {
                string test = it4.key();
                
                if (!test.compare("pdb")) {
                    vector<string> tab = it4.value();
                    list_pdbs2 = tab;

                    //cout << id << " / " << id2 << endl;
                    for (uint k = 0; k < list_pdbs2.size(); k++) {
                        if (count(list_pdbs.begin(), list_pdbs.end(), list_pdbs2[k])) {
                            is_added = false;
                        }
                        //cout << list_pdbs2[k] << endl;
                    }

                } 
                
            }
            //cout << endl;*/
        }
    
       
        /*for(uint ii = 0; ii < list_pfams.size(); ii++) {
            for (uint jj = 0; jj < list_pfams[ii].size(); jj++) {
                cout << "[" << ii << "][" << jj << "]: " << list_pfams[ii][jj] << endl;
            }
        }*/
        if (is_added) {      
            new_id["pdb"] = list_pdbs;     
            new_motif[id] = new_id;
        }
        new_id.clear();
        //cout << "valeur: " << ite << endl;
        /*for (uint i = 0; i < tab_struc.size() ; i++) {
        cout << "tab_struc[" << i << "]: " << tab_struc[i] << endl << endl;
        } */
    }
    outfile << new_motif.dump(4) << endl;
    outfile.close(); 
}

int main()
{
    string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/bibli_test2.json";
    string jsontest = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark_test.json";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final_test.json";
    delete_redundant_pdb(jsonfile, jsontest, out);
    return 0;
}
    
