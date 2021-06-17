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
            string id_pdb = it3.key();
            string sequence2, structure2;
            //cout << "id2: " << id2 << endl;

            //cout << "id: " << id << " / id2: " << id2 << endl;
            for (auto it4 = js2[id_pdb].begin(); it4 != js2[id_pdb].end(); ++it4) {
                string chain = it4.key();

                for (auto it5 = js2[id_pdb][chain].begin(); it5 != js2[id_pdb][chain].end(); ++it5) {
                    string test = it5.key(); 

                    if (!test.compare("sequence")) {
                        sequence2 = it5.value();
                        //cout << sequence2 << endl;
                        if (!sequence.compare(sequence2) && !structure.compare(structure2)) {
                            //cout << id2 << endl;
                            vector<string> tmp;
                            tmp.push_back(id_pdb);
                            cout << "pdb: " << id_pdb << endl;
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

vector<string> find_components(string sequence, string delimiter) {
    vector<string> list;
    string seq = sequence;
    string subseq;
    uint fin = 0;

    while(seq.find(delimiter) != string::npos) {
        fin = seq.find(delimiter);
        
        subseq = seq.substr(0, fin);
        seq = seq.substr(fin + 1);
        list.push_back(subseq); // new component sequence
        //std::cout << "subseq: " << subseq << endl;
    } 
    if (!seq.empty()) {
        list.push_back(seq);
        //std::cout << "subseq: " << seq << endl;
    }
    return list;
}

string is_include(vector<string>& components, string sequence, vector<string>& contacts) {
    
    string seq_contact = "";
    vector<uint> positions;
    uint count = 0;
    uint debut = 0;
    string str = components[0];
    
    uint pos = sequence.find(str, 0);
    debut = pos + components[0].size();

    if (pos == 0) {
        seq_contact += contacts[0];
    } else if (pos <= sequence.size()) {
        string gap = "";
        for (uint i = 0; i < pos; i++) {
            gap += ".";
        }
        seq_contact += gap + contacts[0];
    }
    while(pos <= sequence.size() && count < components.size() - 1)
    {   
        string gap = "";
        debut = pos + components[count].size();    
        count++;
        str = components[count];
        pos = sequence.find(str, pos + components[count-1].size());
        
        for (uint i = debut; i < pos; i++) {
            gap += ".";
        }
        seq_contact += gap + contacts[count];
        
    } 
    if (count == components.size() - 1) {
        string gap = "";
        if (seq_contact.size() != sequence.size()) {  
            for (uint i = 0; i < sequence.size() - seq_contact.size(); i++) {
                gap += ".";
            }
        }
        seq_contact += gap;
        return seq_contact;
    } 
    return std::string();
}

//Concatenate the contact field to the motives of the benchmark (which is obtained from the motives library)
string add_contact(const string& jsonbm, const string& jsonmotifs) {
    std::ifstream lib(jsonbm);
    std::ifstream lib2(jsonmotifs);
    string bm2 = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.json";
    std::ofstream outfile (bm2);
    json new_motif;
    json new_id;
    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        string seq_bm;
        string seq_contact;

        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            string test = it2.key();   
            //cout << "test: " << it2.key() << endl;
            if (!test.compare("seq")) {
                seq_bm = it2.value();
                new_id[test] = it2.value();
            } else {
                new_id[test] = it2.value();
            } 
        }
        //cout << "-------begin---------" << endl;
        
        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();
            vector<string> comp;
            vector<string> strucs;
            vector<string> list_pdbs;
            bool flag = false;

            //cout << "id: " << id << " / id2: " << id2 << endl;
            for (auto it4 = js2[id2].begin(); it4 != js2[id2].end(); ++it4) {
                string test = it4.key();
               
                if (!test.compare("sequence")) {
                    string sequence = it4.value();
                    comp = find_components(sequence, "&");
                    //cout << id << " / " << id2 << endl;
                } else if (!test.compare("contacts")) {
                    string struc2d = it4.value();
                    strucs = find_components(struc2d, "&");
                    
                } else if (!test.compare("pdb")) {
                    vector<string> tab = it4.value();
                    list_pdbs = tab;
                    if (find(list_pdbs.begin(), list_pdbs.end(), id) != list_pdbs.end()) {
                        flag = true;
                    }
                } 
            }
            if (flag) {
                seq_contact = is_include(comp, seq_bm, strucs);
                //cout << "id: " << id << " id2: " << id2 << " seq_contact: " << seq_contact << endl;
                new_id["ctc"] = seq_contact;
            }  
            
        }
    
        new_motif[id] = new_id;
        new_id.clear();

    }
    outfile << new_motif.dump(4) << endl;
    outfile.close(); 
    return bm2;
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
                    outfasta << seq.substr(0,seq.size()) << endl;
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
    //string jsonbeta = path + "modules/ISAURE/Motifs_version_initiale/motifs_beta.json";
    string fasta = path + "fasta/benchmark.fa";
    //string jsondssr = path + "modules/ISAURE/Motifs_version_initiale/dssr2.json";
    string jsonmotifs = path + "modules/ISAURE/Motifs_version_initiale/motifs_06-06-2021.json";
    string jsonbm1 = path + "modules/ISAURE/Motifs_version_initiale/benchmark_16-06-2021.json";

    //string jsonremove = create_benchmark(jsonbeta, fasta, 20);
    //string jsonremovewithpdb = add_pdb(jsonremove, jsondssr);
    string jsonbm2 = add_contact(jsonbm1, jsonmotifs);
    //string jsonlibrary = delete_redundant_pdb(jsonmotifs, jsonremovewithpdb);

    return 0;
}
    
