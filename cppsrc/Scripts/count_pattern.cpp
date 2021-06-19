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

//Return true if the first sequence seq1 is included in the second sequence seq2
//if not return false
int is_contains(string& seq1, string& seq2) {

    uint size1 = seq1.size();
    uint size2 = seq2.size();
    int index = -1;
    if (size1 > size2) {
        //cout << "size1: " << size1 << ", size2: " << size2 << endl;
        return -1;
    }

    /*cout << "seq1: " << seq1 << endl;
    cout << "seq2: " << seq2 << endl;*/
    index = seq2.find(seq1);
    if (index == -1) {
        return -1;
    } else {
        //cout << "index: " << index << endl;
        return index;
    }
    return -1;
    
}

//If we find the sequence and structure of pattern A in pattern B, we have to concatenate the pfam lists of A and B,
//remove the duplicates, assign this new list of pfam lists to A, and assign as occurrence to A the size of this list.
void counting_occurences(const string& jsonfile, const string& jsonoutfile) {
    std::ifstream lib(jsonfile);
    std::ifstream lib2(jsonfile);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;
    string delimiter = "&";

    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
    //the list of pfam lists of the motif we want to count the inclusion in other motif
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        string test;
        uint occurrences = 0;
        int fin;
        string sequence;
        string struc;
        vector<string> composantes;
        vector<string> tab_struc;
        vector<vector<string>> list_pfams;
        vector<vector<string>> list_pfams2;
        vector<vector<string>> union_pfams;
        bool is_change = false;

        //cout << "id: " << id << endl;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            test = it2.key();
            string test = it2.key();   

            if (!test.compare("pfam")) {
                vector<vector<string>> tab = it2.value();
                list_pfams = tab;
                /*set<set<string>>::iterator iit;
                set<string>::iterator iit2;
                for(iit = list_pfams.begin(); iit != list_pfams.end(); iit++) {
                    for (iit2 = iit->begin(); iit2 != iit->end(); ++iit2) {
                        cout << *iit2 << endl;
                    }
                    cout << endl << endl;
                }*/
            } else if (!test.compare("sequence")) {
                //cout << "sequence: " << it2.value() << endl;
                sequence = it2.value();
                new_id[test] = it2.value();

                string subseq;
                while(sequence.find(delimiter) != string::npos) {
                    fin = sequence.find(delimiter);
                    
                    subseq = sequence.substr(0, fin);
                    sequence = sequence.substr(fin + 1);
                    composantes.push_back(subseq); // new component sequence
                    //std::cout << "subseq: " << subseq << endl;
                } 
                if (!sequence.empty()) {
                    composantes.push_back(sequence);
                    //std::cout << "subseq: " << seq << endl;
                }
  
            } else if (!test.compare("struct2d")) {
                //cout << "struct2d: " << it2.value() << endl;
                struc = it2.value();
                new_id[test] = it2.value();
                string subseq;
                while(struc.find(delimiter) != string::npos) {
                    fin = struc.find(delimiter);
                    
                    subseq = struc.substr(0, fin);
                    struc = struc.substr(fin + 1);
                    tab_struc.push_back(subseq); // new component sequence
                    //std::cout << "subseq: " << subseq << endl;
                } 
                if (!struc.empty()) {
                    tab_struc.push_back(struc);
                    //std::cout << "subseq: " << seq << endl;
                }
            }
             else if (!test.compare("occurences") ) {
                occurrences = it2.value();
            } else {
                new_id[test] = it2.value();
            }  
        }
        //cout << "-------begin---------" << endl;
        
        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();
            string sequence2, struc2;
            vector<string> composantes2;
            vector<string> tab_struc2;
            int occurences2;
            int fin;

            //cout << "id: " << id << " / id2: " << id2 << endl;
            for (auto it4 = js[id2].begin(); it4 != js[id2].end(); ++it4) {
                string test = it4.key();
                
                if (id != id2) {       
                    if (!test.compare("pfam")) {
                        vector<vector<string>> tab = it4.value();
                        list_pfams2 = tab;
                        /*for (uint k = 0; k < tab2.size(); k++) {
                            for (uint l = 0; l < tab2[k].size(); l++) {
                                pfams2.insert(tab2[k][l]);
                            }
                            list_pfams2.insert(pfams);
                            pfams2.clear();
                        }*/
            
                        /*set<set<string>>::iterator iit;
                        set<string>::iterator iit2;
                        for(iit = list_pfams.begin(); iit != list_pfams.end(); iit++) {
                            for (iit2 = iit->begin(); iit2 != iit->end(); ++iit2) {
                                cout << *iit2 << endl;
                            }
                            cout << endl << endl;
                        }*/
                    } else if (!test.compare("occurences")) {
                        occurences2 = it4.value();
                        //cout << "occurences2: "<< occurences2 << endl;

                    } else if (!test.compare("sequence")) {
                        sequence2 = it4.value();
                        string subseq;
                        while(sequence2.find(delimiter) != string::npos) {
                            fin = sequence2.find(delimiter);
                            
                            subseq = sequence2.substr(0, fin);
                            sequence2 = sequence2.substr(fin + 1);
                            composantes2.push_back(subseq); // new component sequence
                            //std::cout << "subseq: " << subseq << endl;
                        } 
                        if (!sequence2.empty()) {
                            composantes2.push_back(sequence2);
                            //std::cout << "subseq: " << seq << endl;
                        }

                    } else if (!test.compare("struct2d")) {
                        struc2 = it4.value();
                        string subseq;
                        while(struc2.find(delimiter) != string::npos) {
                            fin = struc2.find(delimiter);
                            
                            subseq = struc2.substr(0, fin);
                            struc2 = struc2.substr(fin + 1);
                            tab_struc2.push_back(subseq); // new component sequence
                            //std::cout << "subseq: " << subseq << endl;
                        } 
                        if (!struc.empty()) {
                            tab_struc2.push_back(struc2);
                            //std::cout << "subseq: " << seq << endl;
                        }
                        
                        uint number = 0;
                        int tab[composantes.size()];
                        for (uint ii = 0; ii < composantes.size(); ii++) {
                            //cout << "tab[" << ii << "]: " << tab[ii] << endl;
                            tab[ii] = 0;
                        }
                        //flag is true if the first component is found or if the k component is indeed placed after the k-1 component
                        //It checks if the found components are in the correct order
                        for (uint k = 0; k < composantes.size() ; k++) {
                            bool flag = false;
                            for (uint l = 0; l < composantes2.size(); l++) {
                                int test1 = is_contains(composantes[k], composantes2[l]);
                                int test2 = is_contains(tab_struc[k], tab_struc2[l]);
                                    if (test1 == test2 && test1 != -1 && test2 != -1) {
                                        if(!flag) {
                                            if (k == 0 || test1 + composantes[k].size() > tab[k-1]) {
                                                tab[k] = test1 + composantes[k].size();
                                                flag = true;
                                            }
                                                
                                        }    
                                    }
                                    //cout << "----end----" << endl;
                                //}
                            }
                            if(flag) {
                               number++;
                            }
                        }
                        
                        // if number equal to the size of the number of component in the motif, it means that the motif is included.
                        //So we add the intersection of the two pfams list to the motif
                        if(number == composantes.size()) {
                            cout << "id: " << id << " / id2: " << id2 << endl;
                            vector<vector<string>> add_pfams;
                            std::set_difference(list_pfams2.begin(), list_pfams2.end(), list_pfams.begin(), list_pfams.end(),
                            std::inserter(add_pfams, add_pfams.begin()));
                            list_pfams.insert(list_pfams.begin(), add_pfams.begin(), add_pfams.end());
                            cout << "size: " << list_pfams.size() << endl;
                            add_pfams.clear();
                            is_change = true;
                        } 
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
        
        new_id["occurences"] = list_pfams.size();
        new_id["pfam"] = list_pfams;
                        
        //cout << "-------ending---------" << endl;
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
    //183
    //cout << "------------------BEGIN-----------------" << endl;
    string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_06-06-2021.json";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_final.json";
    counting_occurences(jsonfile, out);

    //cout << "------------------END-----------------" << endl;
    return 0;
}
    
