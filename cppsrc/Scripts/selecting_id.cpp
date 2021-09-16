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
That script will remove from the library all the pattern that match ONLY with the sequence from which it comes from.
*/

vector<string> get_list_pdb_benchmark(const string& benchmark) {

    ifstream bm(benchmark);
    vector<string> list_pdb;
    if (bm.is_open()) {
        string name;
        string sequence;
        string structure;
        string contacts;

        while (getline(bm, name)) {
            int size = name.size();
            name = name.substr(5,size-8); 
            list_pdb.push_back(name);

            getline(bm, sequence);
            getline(bm, structure);
            getline(bm, contacts);
        }
        bm.close();
    }
    return list_pdb;
}

string trim(string str) {
    int size = str.size();
    str = str.substr(1, size-2);
    return str;
}

bool find_id_pattern(string& pdb_pattern, const string& benchmark) {
    vector<string> l = get_list_pdb_benchmark(benchmark);
    for (string pdb_bm : l) {
        if (!pdb_bm.compare(pdb_pattern)) {
            return true;
        }
    }
    return false;
}

vector<pair<string, string>> find_id(const string& bibli, const string& benchmark) {
    std::ifstream lib(bibli);
    json js = json::parse(lib);

    vector<pair<string, string>> association;
    
    for (auto it = js.begin(); it != js.end(); ++it) {  
        string id = it.key();
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) { 
            string field = it2.key();
            if (!field.compare("pdb")) {
                int n = js[id][field].size();
                for (int i = 0; i < n ; i++) {
                    ostringstream stream;
                    stream << js[id][field][i];
                    string pdb = trim(stream.str());
                    if (find_id_pattern(pdb, benchmark)) {
                        pair<string, string> p;
                        p.first = pdb;
                        p.second = id;
                        association.push_back(p);
                    }
                }
            }
        }
    }
    
    lib.close();
    return association;
}

int main()
{
    string bibli = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final.json";
    string benchmark = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.dbn";

    vector<pair<string, string>> association = find_id(bibli, benchmark);
    /*for (pair<string,string> p : association) {
        cout << "<" << p.first << ", " << p.second << ">" << endl;
    }*/

    return 0;
}