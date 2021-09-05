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
/*
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
}*/
/*
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
}*/

void create_benchmark(const string& jsonmotifs) {
    std::ifstream lib(jsonmotifs);
    string fasta = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/fasta/";
    string list = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.txt";
    string dbn = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.dbn";
    std::ofstream outlist (list);
    std::ofstream outdbn (dbn);
    json js = json::parse(lib);
    uint count = 0;

    for (auto it = js.begin(); it != js.end(); ++it) {    
        string id = it.key();
        string name, seq, contacts, structure;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            string chain = it2.key();
            if (chain.compare("pfams") != 0) {
                string name = id + "_" + chain;
                string filename = fasta + name + ".fa";
                std::ofstream outfasta (filename);
                outfasta << ">test_" << name << endl;
                for (auto it3 = js[id][chain].begin(); it3 != js[id][chain].end(); ++it3) {     
                    string field = it3.key();
                    if (!field.compare("sequence")) {
                        seq = it3.value();
                        outfasta << seq.substr(0,seq.size()) << endl;
                        outfasta.close();

                    } else if (!field.compare("contacts")) {
                        contacts = it3.value();

                    } else if (!field.compare("struct2d")) {
                        structure = it3.value();
                    }
                }
                if(seq.find('&') == string::npos) {
                    outlist << ">test_" << name << endl;
                    outdbn << "test_" << name << "." << endl;
                    outlist << contacts << endl;
                    outdbn << seq << endl;
                    outdbn << structure << endl;
                    outdbn << contacts << endl;
                    outlist << seq << endl;
                    outlist << structure << endl;      
                    count++;       
                }
            }
        }
    }
    cout << count << " sequences en tout" << endl;
    lib.close();
    outlist.close();
    outdbn.close();
}

int main()
{
    string path = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/";
    //string jsonmotifs = path + "modules/ISAURE/Motifs_version_initiale/motifs_06-06-2021.json";
    string jsonbm = path + "modules/ISAURE/Motifs_version_initiale/benchmark_16-07-2021.json";

    
    //string jsonbm2 = add_contact(jsonbm1, jsonmotifs);
    create_benchmark(jsonbm);

    return 0;
}
    
