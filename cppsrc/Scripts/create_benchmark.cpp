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

string delete_redundant_pdb(const string& jsonfile, const string& jsontest) {
    std::ifstream lib(jsonfile);
    std::ifstream lib2(jsontest);
    
    string jsonlibrary = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final.json";

    std::ofstream outfile (jsonlibrary);
    json new_motif;
    json new_id;
    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
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
            } else {
                new_id[test] = it2.value();
            }  
        }       
        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();

            /*if (id == "899")
                cout << "id: " << id << " / id2: " << id2 << endl;*/
            for (auto it4 = js2[id2].begin(); it4 != js2[id2].end(); ++it4) {
                string test = it4.key();
                /*if (id == "899")
                    cout << "id: " << id << " / id2: " << id2 << endl;*/
                if (!test.compare("pdb")) {
                    vector<string> tab = it4.value();
                    list_pdbs2 = tab;
                    for (uint k = 0; k < list_pdbs2.size(); k++) {
                        /*if (id == "142" && id2 == "103") {
                                cout << list_pdbs2[k] << endl;
                                for (uint ii = 0; ii < list_pdbs.size(); ii++) {
                                    cout << list_pdbs[ii] << endl;
                                }
                            }*/
                        if (count(list_pdbs.begin(), list_pdbs.end(), list_pdbs2[k])) {
                            
                            is_added = false;         
                        }
                        //cout << list_pdbs2[k] << endl;
                    }
                } 
                
            }
        }

        if (is_added) {      
            new_id["pdb"] = list_pdbs;     
            new_motif[id] = new_id;
        }
        new_id.clear();
    }
    outfile << new_motif.dump(4) << endl;
    outfile.close(); 
    return jsonlibrary;
}


//Concatenate the motives from jsonmotifs by adding the corresponding pdb from jsondssr
string add_pdb(const string& jsonmotifs, const string& jsondssr) {
    std::ifstream lib(jsonmotifs);
    std::ifstream lib2(jsondssr);
    string jsonwithpdb = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/remove_motif_pdb.json";
    std::ofstream outfile (jsonwithpdb);
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
            //cout << "id2: " << id2 << endl;

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
        }
    
        new_motif[id] = new_id;
        new_id.clear();

    }
    outfile << new_motif.dump(4) << endl;
    outfile.close(); 
    return jsonwithpdb;
}

string find_id_max(const string& jsonmotifs, vector<string>& v) {

std::ifstream lib(jsonmotifs);
json js = json::parse(lib);
size_t max = 0;
string id_max;
bool flag = true;

for (auto it = js.begin(); it != js.end(); ++it) {
    string id = it.key();
    string structure;
    
    if (find(v.begin(), v.end(), id) == v.end()) {
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            string test = it2.key(); 

            if (!test.compare("sequence")) {
                string sequence = it2.value();
                size_t found = sequence.find("&");
                if (found == string::npos && flag) {
                    if (sequence.size() > max) {
                        max = sequence.size();
                        id_max = id;
                    }
                }
            } /*else if (!test.compare("struct2d")) {
                string structure = it2.value();
                cout << "id: " << id << endl;
                size_t found = structure.find("[");
                if (found != string::npos) {
                    flag = false;
                }

            }*/
        }
    }
}   

return id_max;

}

vector<string> vector_id_max(const string& jsonmotifs, vector<string>& v, size_t number) {

std::ifstream lib(jsonmotifs);
json js = json::parse(lib);
size_t count = 0;
string id_max;
vector<string> tab_max;

while (count != number) {
    id_max = find_id_max(jsonmotifs, tab_max);
    tab_max.push_back(id_max);
    count ++;
}

return tab_max;

}

string create_benchmark(const string& jsonmotifs, const string& fasta, size_t number_seq) {
    std::ifstream lib(jsonmotifs);
    
    std::ofstream outfasta (fasta);

    string jsonremove = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/remove_motif.json";
    std::ofstream outmotif (jsonremove);
    json new_motif;
    json new_id;
    json js = json::parse(lib);
    
    vector<string> tab_max = vector_id_max(jsonmotifs, tab_max, number_seq);

    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();

        if (find(tab_max.begin(), tab_max.end(), id) != tab_max.end()) {
            outfasta << ">test_" << id << endl;
            for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
                string test = it2.key();
                if (!test.compare("sequence")) {
                    string seq = it2.value();
                    outfasta << seq.substr(0,seq.size()-1) << endl;
                }
                new_id[test] = it2.value();   
            }
            new_motif[id] = new_id;
            new_id.clear();
        }
    }
    outmotif << new_motif.dump(4) << endl;
    outmotif.close(); 
    return jsonremove;
}

int main()
{
    string path = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/";
    string jsonbeta = path + "modules/ISAURE/Motifs_version_initiale/motifs_beta.json";
    string fasta = path + "fasta/benchmark.fasta";
    string jsondssr = path + "modules/ISAURE/Motifs_version_initiale/dssr2.json";
    string jsonmotifs= path + "modules/ISAURE/Motifs_version_initiale/motifs_06-06-2021.json";

    string jsonremove = create_benchmark(jsonbeta, fasta, 20);
    string jsonremovewithpdb = add_pdb(jsonremove, jsondssr);
    string jsonlibrary = delete_redundant_pdb(jsonmotifs, jsonremovewithpdb);

    return 0;
}
    
