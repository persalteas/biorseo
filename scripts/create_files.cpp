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
Create a .fasta file for each of the sequence inside the benchmark in json format.
Also create a .dbn and .txt file that list the name, sequence, 2d structure and contacts for all sequence in the benchmark file.
Those files are useful for the Isaure_benchmark.py script.
*/
void create_files(const string& jsonmotifs) {
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
    string jsonbm = path + "modules/ISAURE/benchmark_16-07-2021.json"; 
    create_files(jsonbm);

    return 0;
}
    
