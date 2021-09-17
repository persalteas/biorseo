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

    fstream bm(benchmark);
    vector<string> list_pdb;
    if (bm.is_open()) {
        string name;
        string sequence;
        string structure;
        string contacts;

        while (getline(bm, name)) {
            int size = name.size();
            name = name.substr(5,size-6); 
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

string find_id_pattern(string& pdb_pattern, const string& benchmark) {
    vector<string> l = get_list_pdb_benchmark(benchmark);
    for (string pdb_bm : l) {
        int size = pdb_bm.size();
        string cmp = pdb_bm.substr(0, size-2);
        if (!cmp.compare(pdb_pattern)) {
            return pdb_bm;
        }
    }
    return string();
}

vector<pair<string, string>> find_id(const string& bibli, const string& benchmark) {
    ifstream lib(bibli);
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
                    string pdb_complete = find_id_pattern(pdb, benchmark);
                    if (!(pdb_complete.empty())) {
                        pair<string, string> p;
                        p.first = pdb_complete;
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

bool does_it_match(const string& result, const string& id_motif) {
    ifstream f_res(result);
    if (f_res.is_open()) {
        string name;
        string seq;
        string struc;
        string contacts;

        getline(f_res, name);
        getline(f_res, seq);
        while (getline(f_res, struc)) {
            string motif_json = "JSON" + id_motif + " +";
            if(struc.find(motif_json, 0) != string::npos) {
                return true;
            }
            motif_json = "JSON" + id_motif + "\n";
            if(struc.find(motif_json, 0) != string::npos) {
                return true;
            }
            getline(f_res,contacts);
        }
        f_res.close();
    }
    return false;
}

vector<string> select_not_motif(const string& bibli, const string& benchmark) {
    vector<string> selection;
    vector<pair<string, string>> association = find_id(bibli, benchmark);
    vector<string> list_bm = get_list_pdb_benchmark(benchmark);

    string path_begin = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/results/test_";
    string path_MFE_F = ".json_pmF_MEA";

    for (pair<string, string> p : association) {
        string id_motif = p.second;
        selection.push_back(id_motif);
    }
    for (pair<string, string> p : association) {
        cout << p.first << ", " << p.second << endl;
    }
    cout << "size: " << association.size() << endl;
   
    for (string pdb : list_bm) {
        string path_result = path_begin + pdb + path_MFE_F;
        for (pair<string,string> pair : association) {
            if (pair.first.substr(0, pair.first.size()-2).compare(pdb.substr(0, pdb.size()-2)) != 0) {
                bool test = does_it_match(path_result, pair.second);

                if (test) {
                    //if (!(pair.second.compare("954"))) { cout << "p1: " << pair.first << "pdb: " << pdb << endl;}
                    auto position = find(selection.begin(), selection.end(), pair.second);
                    if (position != selection.end()) {
                        int index = position - selection.begin();
                        selection.erase(selection.begin() + index);
                    }
                }
            }
        }
    }
    sort(selection.begin(), selection.end() );
    selection.erase(unique(selection.begin(), selection.end() ), selection.end() );

    cout << "size: " << selection.size() << endl;

    return selection;
}

int main()
{
    string bibli = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final.json";
    string benchmark = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.dbn";

    /*vector<pair<string, string>> association = find_id(bibli, benchmark);
    for (pair<string,string> p : association) {
        cout << "<" << p.first << ", " << p.second << ">" << endl;
    }*/

    vector<string> selection = select_not_motif(bibli, benchmark);
    for (string str : selection) {
        cout << str << ", ";
    }
    cout << endl;

    /*string result = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/results/test_1U6P_B.json_pmF_MEA";
    bool test = does_it_match(result, "150");
    cout << "test : " << test << endl;*/

    return 0;
}