#include <iostream>
#include <sstream>
#include <fstream>
#include "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/cppsrc/json.hpp"
#include <typeinfo>
#include <set>
#include <algorithm>

using namespace std;
using json = nlohmann::json;

//Create a new file base on the one in argument that will contains one set of pfams for each pattern
string pfams_union (const string& jsonfile) {
    std::ifstream lib(jsonfile);
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/test.json";
    std::ofstream outfile (out);

    json js = json::parse(lib);
    json new_motif;
    json new_id;

    std::string keys[5] = {"occurences", "pdb", "pfam", "sequence", "struct2d"};
    for (auto it = js.begin(); it != js.end(); ++it) {
        int j = 0;
        string id = it.key();
        //std::vector<std::string> pfams = js["pfam"];

         for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {
            string test = it2.key();
            set<string> pfams;
            //json js2 = it.value();
            if (!test.compare(keys[2])) {
                vector<vector<string>> tab = it2.value();
                for (uint k = 0; k < tab.size(); k++) {
                    for (uint l = 0; l < tab[k].size(); l++) {
                        pfams.insert(tab[k][l]);
                        //cout << "names[" << k << "][" << l << "]: " << tab[k][l] << endl;
                    }
                }
                new_id[test] = pfams;
                /*cout << endl;
                set<string>::iterator ii ;
                cout << "The element of set of size " << pfams.size() << " are : \n";
                for (ii = pfams.begin() ; ii != pfams.end() ; ii++ ) 
                {
                    cout << *ii<<" ";
                }
                cout << endl;*/
            } else {
                new_id[test] = it2.value();
            }       
        }
        new_motif[id] = new_id;
        new_id.clear();
    }
    //new_js.push_back(new_id);
    outfile << new_motif.dump(4) << endl;
    outfile.close();
    return out;
}

bool is_contains_set(set<string>& s1, set<string>& s2) {
    //cout << "-----begin------" << endl;
    set<string>::iterator subset;
    set<string>::iterator it ;

    uint size1 = s1.size();
    uint size2 = s2.size();
    
    if (size1 > size2) {
        //cout << "size1: " << size1 << ", size2: " << size2 << endl;
        return false;
    }
    for(string s: s1) {
        //cout << "s1: " << s << endl;
        if(s2.count(s) == 0) {
            //cout << "count: " << s2.count(s) << endl;
            return false;
        }
    }
    //cout << "-----end------" << endl;
    return true;
}

void counting_occurences(const string& jsonfile) {
    std::ifstream lib(jsonfile);
    std::ifstream lib2(jsonfile);
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/motifs.json";
    std::ofstream outfile (out);
    set<string> pfams;
    set<string> pfams2;

    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        uint count = 0;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {
            string test = it2.key();
            
            if (!test.compare("pfam")) {
                vector<string> tmp = it2.value();
                //cout << "pfams: " << endl;
                for (uint i = 0; i < tmp.size(); i++) {
                    pfams.insert(tmp[i]);
                    //cout << tmp[i] << endl;
                }  
            }
        }

        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();
            for (auto it4 = js[id2].begin(); it4 != js[id2].end(); ++it4) {
                string test = it4.key();
                
                if (id != id2 && !test.compare("pfam")) {
                    vector<string> tmp2 = it4.value();
                    //cout << "pfams2: " << endl;
                    for (uint i = 0; i < tmp2.size(); i++) {
                        pfams2.insert(tmp2[i]);
                        //cout << tmp2[i] << endl;
                    }  
                    if (is_contains_set(pfams, pfams2)) {
                        count++;
                        //cout << id << "/" << id2 << ": " << count << endl;
                    }
                }

            }
            //cout << endl;
            pfams2.clear();
        }
        pfams.clear();
        //cout << "valeur: " << ite << endl;    
    }

    outfile << "coucou" << endl;
    outfile.close();
}

int main()
{
    cout << "------------------BEGIN-----------------" << endl;
    string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_28-05-2021.json";
    string tmpfile = pfams_union(jsonfile);
    //cout << tmpfile << endl;
    counting_occurences(tmpfile);

    cout << "------------------END-----------------" << endl;
}