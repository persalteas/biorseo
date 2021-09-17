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

struct data { 
    string pdb;
    string seq_pdb;
    string id;
    string cmp;
};
typedef struct data data;


vector<data> get_list_pdb_benchmark(const string& benchmark) {

    fstream bm(benchmark);
    vector<data> list_pdb_seq;
    if (bm.is_open()) {
        string name;
        string sequence;
        string structure;
        string contacts;

        while (getline(bm, name)) {
            data d;
            int size = name.size();
            name = name.substr(5,size-6); 
            getline(bm, sequence);
            d.pdb = name;
            d.seq_pdb = sequence;
            list_pdb_seq.push_back(d);

            getline(bm, structure);
            getline(bm, contacts);
        }
        bm.close();
    }
    return list_pdb_seq;
}

string trim(string str) {
    int size = str.size();
    str = str.substr(1, size-2);
    return str;
}

data find_id_pattern(string& pdb_pattern, const string& benchmark) {
    vector<data> l = get_list_pdb_benchmark(benchmark);
    int size = l.size();

    for (data d : l) {
        string cmp = d.pdb;
        cmp = cmp.substr(0, d.pdb.size()-2);
        if (!cmp.compare(pdb_pattern)) {
            return d;
        }
    }
    return data();
}

vector<data> find_id(const string& bibli, const string& benchmark) {
    ifstream lib(bibli);
    json js = json::parse(lib);

    //nam seq_bm et id seq_id
    vector<data> association;
    
    for (auto it = js.begin(); it != js.end(); ++it) {  
        string id = it.key();
        data d;

        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) { 
            string field = it2.key();
            string seq;
            if (!field.compare("pdb")) {
                int n = js[id][field].size();
                for (int i = 0; i < n ; i++) {
                    ostringstream stream;
                    stream << js[id][field][i];
                    string pdb = trim(stream.str());
                    
                    d = find_id_pattern(pdb, benchmark);
                }
            }

            if (!field.compare("sequence")) {
                seq = it2.value();

                if (!(d.pdb.empty())) {                    
                    d.id = id;
                    d.cmp = seq;
                    association.push_back(d);
                }
            }
        }
    }
    lib.close();
    cout << association.size() << endl;
    return association;
}

bool does_it_match(const string& seq, const string& seq_motif) {
    size_t found = seq_motif.find("&");
    size_t size = seq_motif.size();
    vector<string> list_cmp;
    if (found != std::string::npos) {
        int count = 1;
        
        string cmp = seq_motif.substr(0, found);
        list_cmp.push_back(cmp);
        while(found != std::string::npos) {
            size_t begin = found;
            found = seq_motif.find("&", found + 1);
            cmp = seq_motif.substr(begin+1, found-begin-1);
            list_cmp.push_back(cmp);
            count++;
        }

        found = seq.find(list_cmp[0]);
        int count2 = 1;
        while((found != std::string::npos) && (count2 < count)) {
            size_t begin = found;
            found = seq.find(list_cmp[count2], found + 1);
            count2++;
        }

        if(count == count2) {
            return true;
        }

    } else {
        found = seq.find(seq_motif);
        if (found != std::string::npos) {
            return true;
        }
    }
    return false;
}

vector<string> select_not_motif(const string& bibli, const string& benchmark) {
    vector<string> selection;
    vector<data> association = find_id(bibli, benchmark);

    for (data d : association) {
        selection.push_back(d.id);
    }

    for (data d : association) {
        for (data d2 : association) {
            string seq = d.seq_pdb;
            string seq2 = d2.cmp;
            bool test = false;

            if(d.pdb.substr(0, d.pdb.size()-2) != d2.pdb.substr(0, d2.pdb.size()-2)) {
                test = does_it_match(seq, seq2);
                if (test) {
                    cout << "pdb: " << d.pdb << " vs " << d2.pdb << " " << d2.cmp << " " << d2.id << endl;
                    auto position = find(selection.begin(), selection.end(), d.id);
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

    /*vector<data> v = get_list_pdb_benchmark(benchmark);
    for (data d : v) {
        cout << d.pdb << ", " << d.seq_pdb << endl;
    }*/

    /*string name = "1U6P_B";
    data d = find_id_pattern(name, benchmark);
    cout << "name: " << d.pdb << ", seq: " << d.seq_pdb << endl;*/

    /*vector<data> association = find_id(bibli, benchmark);
    for (data d : association) {
        cout << "<" << d.pdb << ", " << d.seq_pdb << ">, " << "<" << d.id << ", " << d.cmp << ">" << endl;
    }*/

    /*string seq = "UGCGCUUGGCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUU";
    string seq_motif = "UGCGCUUGGCGUUUUAGAGC&GCAAGUUAAAAUAAGGCUAGUCCGUUAUCAA&UGGCACCGAGUCG&U";
    bool test = does_it_match(seq, seq_motif);
    cout << test << endl;*/

    vector<string> selection = select_not_motif(bibli, benchmark);
    for (string str : selection) {
        cout << str << ", ";
    }
    cout << endl;

    return 0;
}