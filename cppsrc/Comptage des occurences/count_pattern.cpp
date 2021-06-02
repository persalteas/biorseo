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

//Create a new file base on the one in argument that will contains one set of pfams for each pattern
string pfams_union (const string& jsonfile) {
    std::ifstream lib(jsonfile);
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/test.json";
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
                //cout << "id: " << id << endl;
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

void counting_occurences2(const string& jsonfile, const string& jsonoutfile) {
    std::ifstream lib(jsonfile);
    std::ifstream lib2(jsonfile);
    
    std::ofstream outfile (jsonoutfile);
    set<string> pfams;
    set<string> pfams2;
    json new_motif;
    json new_id;

    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        string test;
        uint count = 0;
        uint occurrences = 0;
        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {
            test = it2.key();
            
            //If we desire to count occurences based on pdb rather than pfam just replace "pfam" by "pdb"
            if (!test.compare("pdb")) {
                vector<string> tmp = it2.value();
                //cout << "pfams: " << endl;
                for (uint i = 0; i < tmp.size(); i++) {
                    pfams.insert(tmp[i]);
                    new_id[test] = it2.value();
                    //cout << tmp[i] << endl;
                }
                    
            } else if (!test.compare("occurences")) {
                occurrences = it2.value();
            } else {
                new_id[test] = it2.value();
            }  
        }
        //cout << "-------begin---------" << endl;
        
        for (auto it3 = js2.begin(); it3 != js2.end(); ++it3) {
            string id2 = it3.key();
            for (auto it4 = js[id2].begin(); it4 != js[id2].end(); ++it4) {
                string test = it4.key();
                
                //If we desire to count occurences based on pdb rather than pfam just replace "pfam" by "pdb"
                if (id != id2 && !test.compare("pdb")) {
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
        new_id["occurences"] = count + occurrences;
        
        //cout << "-------ending---------" << endl;
        pfams.clear();
        new_motif[id] = new_id;
        new_id.clear();
        //cout << "valeur: " << ite << endl;    
    }

    outfile << new_motif.dump(4) << endl;
    outfile.close();
}

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

void counting_occurences(const string& jsonfile, const string& jsonoutfile) {
    std::ifstream lib(jsonfile);
    std::ifstream lib2(jsonfile);
    
    std::ofstream outfile (jsonoutfile);
    json new_motif;
    json new_id;
    string delimiter = "&";

    json js = json::parse(lib);
    json js2 = json::parse(lib2);
    
    for (auto it = js.begin(); it != js.end(); ++it) {
        string id = it.key();
        string test;
        uint count = 0;
        uint occurrences = 0;
        int fin;
        string sequence;
        string struc;
        vector<string> composantes;
        vector<string> tab_struc;

        for (auto it2 = js[id].begin(); it2 != js[id].end(); ++it2) {      
            test = it2.key();
            
            if (!test.compare("sequence")) {
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
             else if (!test.compare("occurences")) {
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
            int fin;
            //cout << "id: " << id << " / id2: " << id2 << endl;
            for (auto it4 = js[id2].begin(); it4 != js[id2].end(); ++it4) {
                string test = it4.key();
                
                if (id != id2) {
                    if (!test.compare("sequence")) {
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
                        /*if(id == "1" && id2 == "100006")
                            for (uint j = 0; j < tab_struc2.size() ; j++) {
                                cout << "tab_struc2[" << j << "]: " << tab_struc2[j] << endl << endl;
                                cout << "composantes2[" << j << "]: " << composantes2[j] << endl << endl;
                            }*/
                        
                        uint number = 0;
                        int tab[composantes.size()];
                        for (uint k = 0; k < composantes.size() ; k++) {
                            bool flag = false;
                            for (uint l = 0; l < composantes2.size(); l++) {

                                int test1 = is_contains(composantes[k], composantes2[l]);
                                int test2 = is_contains(tab_struc[k], tab_struc2[l]);
                                /*if (id == "1" && id2 == "100006") {
                                    
                                    cout << endl;
                                }*/
                                /*if (id == "10" && id2 == "65") {
                                    cout << "----begin----" << endl;
                                    cout << "id: " << id << " id2: " << id2 << endl;
                                    cout << "[" << k << "][" << l << "]: " << endl;
                                    cout << "test1: " << test1 << " | " << "test2: " << test2 << endl;
                                    cout << "seq1: " << composantes[k] << " vs " << "seq2: " << composantes2[l] << endl;
                                    cout << "struc1: " << tab_struc[k] << " vs " << "struc2: " <<tab_struc2[l] << endl << endl;*/
                                if (test1 == test2 && test1 != -1 && test2 != -1) {
                                    //cout << "flag: " << flag << endl;
                                    if(!flag) {
                                        if (k == 0 || test1 > tab[k-1])
                                            tab[k] = test1;
                                            flag = true;
                                    }    
                                }
                                    //cout << "----end----" << endl;
                                //}
                            }
                            if(flag) {
                               number++;
                            }
                        }
                        if(number == composantes.size()) {
                            count++;
                            cout << id << " vs " << id2 << endl;
                        }
                        /*int test1 = is_contains(sequence, sequence2);
                        int test2 = is_contains(struc, struc2);
                        if (test1 == test2 && test1 != -1 && test2 != -1) {
                            cout << is_contains(sequence, sequence2) << " / " << is_contains(struc, struc2)<< endl;
                            count++;
                            //cout << id << "/" << id2 << ": " << count << endl;
                        }*/
                    }
                }
            }
            //cout << endl;
            /*for (uint j = 0; j < tab_struc2.size() ; j++) {
                cout << "tab_struc2[" << j << "]: " << tab_struc2[j] << endl << endl;
            }*/  
        }
        new_id["occurences"] = count + occurrences;
        
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
    cout << "------------------BEGIN-----------------" << endl;
    string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/motifs_28-05-2021.json";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final.json";
    string tmpfile = pfams_union(jsonfile);
    counting_occurences(tmpfile, out);

    /*string jsonfile = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/bibli_test.json";
    string out = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version/motifs_final_test.json";
    string tmpfile = pfams_union(jsonfile);
    counting_occurences(tmpfile, out);*/

    /*if (std::remove(tmpfile.c_str()) != 0)
		perror("File deletion failed \n");
	else
		cout << "File deleted successfully" << endl;*/

    /*string s1 = "(..)";
    string s2 = "((..))";
    bool test = is_contains(s1, s2);
    cout << "test: " << test << endl;*/ 

    cout << "------------------END-----------------" << endl;
    return 0;

    
}