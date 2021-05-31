#include <iostream>
#include <sstream>
#include <fstream>
#include "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/cppsrc/json.hpp"
#include <typeinfo>
#include <set>

using namespace std;
using json = nlohmann::json;

int main()
{
    cout << "------------------BEGIN-----------------" << endl;
    string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_28-05-2021.json";
    std::ifstream lib(jsonfile);
    std::ofstream outfile ("/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/test.json");

    json js = json::parse(lib);
    json new_motif;
    json new_id;

    std::string keys[5] = {"occurences", "pdb", "pfam", "sequence", "struct2d"};
    for (auto it = js.begin(); it != js.end(); ++it) {
        int j = 0;
        string ite = it.key();
        //std::vector<std::string> pfams = js["pfam"];

         for (auto it2 = js[ite].begin(); it2 != js[ite].end(); ++it2) {
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
        new_motif[ite] = new_id;
        new_id.clear();
    }
    //new_js.push_back(new_id);
    outfile << new_motif.dump(4) << endl;
    
    outfile.close();
    cout << "------------------END-----------------" << endl;
}