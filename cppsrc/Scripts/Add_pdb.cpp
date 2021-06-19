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

//Concatenate the motives from jsonmotifs by adding the corresponding pdb from jsondssr
void add_pdb(const string& jsonmotifs, const string& jsondssr, const string& jsonoutfile) {
    std::ifstream lib(jsonmotifs);
    std::ifstream lib2(jsondssr);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;
    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();

        string sequence, structure;
        vector<string> list_pdbs;
        vector<string> list_pdbs2;
        bool is_added = true;

        //cout << "id: " << id << endl;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            string test = it2.key();   

            if (!test.compare("sequence")) {
                sequence = it2.value();
                new_id[test] = it2.value();

            } else if (!test.compare("struct2d")) {
                structure = it2.value();
                new_id[test] = it2.value();
            
            } else {
                new_id[test] = it2.value();
            }  
        }
        //cout << "-------begin---------" << endl;
        
        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();
            string sequence2, structure2;

            //cout << "id: " << id << " / id2: " << id2 << endl;
            for (auto it4 = js2[id2].begin(); it4 != js2[id2].end(); ++it4) {
                string chain = it4.key();

                for (auto it5 = js2[id2][chain].begin(); it5 != js2[id2][chain].end(); ++it5) {
                    string test = it5.key(); 

                    if (!test.compare("sequence")) {
                        sequence2 = it5.value();
                        //cout << sequence2 << endl;
                        if (!sequence.compare(sequence2) && !structure.compare(structure2)) {
                            //cout << id2 << endl;
                            vector<string> tmp;
                            tmp.push_back(id2);
                            new_id["pdb"] = tmp;
                        }

                    } else if (!test.compare("2D      ")) {
                        structure2 = it5.value();
                        //cout << structure2 << endl;
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
        new_motif[id] = new_id;
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
    string jsonmotifs = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_beta.json";
    string jsondssr = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/dssr2.json";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_fusion_beta.json";
    add_pdb(jsonmotifs, jsondssr, out);
    return 0;
}
    
