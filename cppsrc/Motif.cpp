#include "Motif.h"
#include "Pool.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <regex>
#include <sstream>
#include <thread>
#include <json.hpp>

using namespace boost::filesystem;
using namespace std;
using json = nlohmann::json;

struct recursive_directory_range {
    typedef recursive_directory_iterator iterator;
    recursive_directory_range(path p) : p_(p) {}

    iterator begin() { return recursive_directory_iterator(p_); }
    iterator end() { return recursive_directory_iterator(); }

    path p_;
};

Motif::Motif(void) {}


Motif::Motif(const vector<Component>& v, string PDB) : comp(v), PDBID(PDB)
{
    is_model_ = false;
    reversed_ = false;
    source_   = RNA3DMOTIF;
}

Motif::Motif(const vector<Component>& v, uint id) : comp(v)
{
    carnaval_id = to_string(id);
    is_model_ = false;
    reversed_ = false;
    source_   = CONTACTS;
}

Motif::Motif(string csv_line)
{
    vector<string> tokens;
    split(tokens, csv_line, boost::is_any_of(","));

    if (csv_line.find(string("True")) != std::string::npos or csv_line.find(string("False")) != std::string::npos) // This has been created by jar3d
    {
        atlas_id = tokens[0];
        score_   = stoi(tokens[2]);
        comp.push_back(Component(make_pair<int, int>(stoi(tokens[3]), stoi(tokens[4]))));
        if (tokens[5] != "-") comp.push_back(Component(make_pair<int, int>(stoi(tokens[5]), stoi(tokens[6]))));
        reversed_ = (tokens[1] == "True");
        is_model_ = true;
        PDBID     = "";
        source_   = RNAMOTIFATLAS;
    }

    else // this has been created by BayesPairing
    {
        score_ = stoi(tokens[1]);

        // identify source:
        if (tokens[0].find(string("rna3dmotif")) == std::string::npos)
        {
            is_model_ = true;
            PDBID     = "";
            source_   = RNAMOTIFATLAS;
            atlas_id  = tokens[0];
        }

        else
        {
            is_model_ = false;
            PDBID     = tokens[0];
            source_   = RNA3DMOTIF;
            atlas_id  = "";
        }

        uint i = 2;
        //while (i < tokens.size())
        while (i < tokens.size()-1)
        {
            if (stoi(tokens[i]) < stoi(tokens[i + 1]))
                comp.push_back(Component(make_pair<int, int>(stoi(tokens[i]), stoi(tokens[i + 1]))));
            else
                i = i + 0;
            i += 2;
        }
    }
}

Motif::Motif(const vector<Component>& v, path rinfile, uint id, bool reversed) : comp(v), reversed_(reversed)
{
    // Loads a motif from the RIN file of Carnaval
    carnaval_id = to_string(id);
    source_     = CARNAVAL;
    is_model_     = false;

    std::ifstream file(rinfile);

    if (file.is_open())
    {
        string line;

        getline(file,line); //skip the header_link line
        getline(file,line); //get the links line
        string link_str;
        size_t index = 0;
        string nt_str;
        size_t sub_index = 0;

        while (line != "")
        {
            Link link;

            //link.nts
            index         = line.find(";");
            link_str     = line.substr(0, index);
            line.erase(0, index+1);

            sub_index     = link_str.find(",");
            nt_str         = link_str.substr(0, sub_index);
            link_str.erase(0, sub_index+1);
            link.nts.first = stoi(nt_str);

            sub_index     = link_str.find(",");
            nt_str         = link_str.substr(0, sub_index);
            link_str.erase(0, sub_index+1);
            link.nts.second = stoi(nt_str);

            //link.long_range
            link.long_range = (link_str == "True");

            links_.push_back(link);
        }


        getline(file,line); //skip the header_comp line
        string pos_str, sub_pos_str, k_str, seq;

        while ( std::getline(file,line) )
        {
            if (line == "\n") break; //skip last line (empty)

            Component c(0,0);

            //c.pos
            index         = line.find(";");
            pos_str     = line.substr(0, index);
            line.erase(0, index+1);

            sub_index     = pos_str.find(",");
            sub_pos_str = pos_str.substr(0, sub_index);
            pos_str.erase(0, sub_index+1);
            c.pos.first = stoi(sub_pos_str);
            c.pos.second = stoi(pos_str);

            //c.k
            index     = line.find(";");
            k_str     = line.substr(0, index);
            line.erase(0, index+1);
            c.k     = stoi(k_str);

            //c.seq_
            c.seq_ = line;

            comp.push_back(c);
        }
    }

    else cout << "\t> RIN file not found : " << rinfile << endl;
}

string Motif::pos_string(void) const
{
    stringstream s;
    s << atlas_id << " ( ";
    for (auto c : comp) s << c.pos.first << '-' << c.pos.second << ' ';
    s << ')';
    return s.str();
}

string Motif::get_identifier(void) const
{
    switch (source_) {
    case RNAMOTIFATLAS: return atlas_id; break;
    case CARNAVAL: return string("RIN") + carnaval_id; break;
    case CONTACTS: return string("JSON") + contacts_id; break;
    default: return PDBID;
    }
}

char Motif::is_valid_DESC(const string& descfile)
{
    // /!\ returns 0 if no errors

    std::ifstream  motif;
    string         line;
    vector<string> bases;
    char           c    = 'a';
    char*          prev = &c;

    motif = std::ifstream(descfile);
    getline(motif, line);    // ignore "id: number"
    getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    for (vector<string>::iterator b = (bases.begin() + 1); b != (bases.end() - 1); b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = stoi(b->substr(0, b->find('_')));

        if (string(1, nt).find_first_not_of("ACGU") != string::npos) return nt;
        if (pos <= 0) return '-';
    }

    while (getline(motif, line)) {
        uint   slash = line.find('/');
        string interaction(line.substr(slash - 1, 3));

        string b1   = line.substr(line.find('(') + 1, 8);    // The first base (position_nucleotide)
        string temp = line.substr(slash + 1);
        string b2   = temp.substr(temp.find('(') + 1, 8);    // The second base (position_nucleotide)
        b1.erase(remove(b1.begin(), b1.end(), ' '), b1.end());
        b2.erase(remove(b2.begin(), b2.end(), ' '), b2.end());
        int p1 = stoi(b1.substr(0, b1.find('_')));
        int p2 = stoi(b2.substr(0, b2.find('_')));

        if ((p2 - p1 != 1) and !interaction.compare("C/C")) return 'b';
        if ((p2 - p1 < 4) and (!interaction.compare("+/+") or !interaction.compare("-/-"))) return 'l';
    }
    return 0;
}

char Motif::is_valid_RIN(const string& rinfile) 
{
    string line;
    uint   n_basepairs, index, motif_length = 0;

    std::ifstream motif = std::ifstream(rinfile);
    getline(motif, line); //skip the header_link line
    getline(motif, line); //get the links line
    n_basepairs = count_if(line.begin(), line.end(), [](char c){ return (c==';'); });
    getline(motif, line); //skip the header_comp line
    while (getline(motif, line))
    {
        // lines are formatted like:
        // pos;k;seq
        // 0,1;2;GU
        if (line == "\n") break; //skip last line (empty)
        size_t start = line.find(";") + 1;
        index = line.find(";", start); // find the second ';'
        motif_length += stoi(line.substr(line.find(";")+1, index));
    }

    if (motif_length < 5) return 'l';

    if (!n_basepairs) return 'x';
    
    return (char) 0;
}

//temporaire---------------------------------------------------

bool checkSecondaryStructure(string struc)
{ 
    stack<uint> parentheses;
    stack<uint> crochets;
    stack<uint> accolades;
    stack<uint> chevrons;
    for (uint i = 0; i < struc.length(); i++)
    {

        if (struc[i] != '(' && struc[i] != ')' 
        && struc[i] != '.' && struc[i] != '&'
        && struc[i] != '[' && struc[i] != ']'
        && struc[i] != '{' && struc[i] != '}'
        && struc[i] != '<' && struc[i] != '>') {
            return false;
        } else {
            for (uint i = 0; i < struc.size(); i++) {
                if (struc[i] == '(') {
                    parentheses.push(i);

                } else if (struc[i] == ')') {
                    if (!parentheses.empty())
                        parentheses.pop();
                    else return false;

                } else if (struc[i] == '[') {
                    crochets.push(i);

                } else if (struc[i] == ']') {
                    if (!crochets.empty())
                        crochets.pop();
                    else return false;

                } else if (struc[i] == '{') {
                    accolades.push(i);

                } else if (struc[i] == '}') {
                    if (!accolades.empty())
                        accolades.pop();
                    else return false;

                } else if (struc[i] == '<') {
                    chevrons.push(i);

                } else if (struc[i] == '>') {
                    if (!chevrons.empty())
                        chevrons.pop();
                    else return false;
                } 
            }
        }
    }
    return (parentheses.empty() && crochets.empty() && accolades.empty() && chevrons.empty());
}


//--------------------------------------------------------------
vector<pair<uint,char>> Motif::is_valid_JSON(const string& jsonfile)
{
    // /!\ returns 0 if no errors
    //cout << "---begin----" << endl;
    std::ifstream  motif;
    motif = std::ifstream(jsonfile);
    json js = json::parse(motif);
    vector<pair<uint,char>> errors_id;
    vector<string> components;
    uint fin = 0;

    std::string keys[5] = {"occurences", "pdb", "pfam", "sequence", "struct2d"};
    for (auto i = js.begin(); i != js.end(); ++i) {
        int j = 0;
        string id = i.key();
        //cout << id << ": " << endl;
        for (auto it = js[id].begin(); it != js[id].end(); ++it) {
            string test = it.key();
            //std::cout << "test: " << test << endl;
            if (test.compare(keys[j])){ 
                //std::cout << "error header : keys[" << j << "]: " << keys[j] << " vs test: " << test << endl;
                errors_id.push_back(make_pair(stoi(id), 'd')); 
                //return 'd'; 
            } else if(!test.compare(keys[4])) {
                //std::cout << "struct2d: " << it.value() << endl;
                string ss = it.value();
                if (ss.empty()) {
                    //std::cout << "error empty" <<endl;
                    errors_id.push_back(make_pair(stoi(id), 'f'));
                    //return 'f';
                } else if (!checkSecondaryStructure(ss)) {
                    //std::cout << "error bracket" <<endl;
                    errors_id.push_back(make_pair(stoi(id), 'n'));
                    //return 'n';
                }
            } else if (!test.compare(keys[3])) {
                //std::cout << "sequence: " << it.value() << "\n";
                /*if (id == "484")
                    cout << "seq: " << it.value() << endl;*/
                string seq = it.value();
                if (seq.empty()) {
                    //std::cout << "error empty 2" <<endl;
                    errors_id.push_back(make_pair(stoi(id), 'e'));
                    //return 'l';
                } else if (seq.size() == 1) {
                    //std::cout << "error too short" << endl;
                    errors_id.push_back(make_pair(stoi(id), 'l'));
                } else {
                string subseq;
                
                    while(seq.find('&') != string::npos) {
                        fin = seq.find('&');  
                        subseq = seq.substr(0, fin);
                        seq = seq.substr(fin + 1);
                        if (subseq.size() >= 2) {
                            components.push_back(subseq); 
                            //std::cout << "subseq: " << subseq << endl;
                        } else {
                            errors_id.push_back(make_pair(stoi(id), 'k'));
                            //std::cout << "error too short1" << endl;
                        }
                    } 
                    if (seq.size() >= 2) {
                        components.push_back(seq);
                        //std::cout << "subseq: " << seq << endl;
                    } else {
                        errors_id.push_back(make_pair(stoi(id), 'k'));
                        //std::cout << "error too short2" << endl;
                    }
                }
            }
            j++;
            //cout << "test fin" << endl << endl;
        }
    //std::cout << "no error!\n" << endl;
    }
    return errors_id;
    //cout << "---end----" << endl;
}
/*char Motif::is_valid_JSON(const string& jsonfile)
{
    // /!\ returns 0 if no errors

    std::ifstream  motif;
    motif = std::ifstream(jsonfile);
    json js = json::parse(motif);

    std::string keys[3] = {"occurences", "sequence", "struct2d"};
    
    for (auto i = js.begin(); i != js.end(); ++i) {
        int j = 0;
        string ite = i.key();
        cout << ite << ": " << endl;
        for (auto it = js[ite].begin(); it != js[ite].end(); ++it) {
            string test = it.key();
            if (test.compare(keys[j])){ 
                return 'd'; 
            }

            if(!test.compare(keys[2])) {
                std::cout << "struct2d: " << it.value() << endl;
                string ss = it.value();
                if (ss.empty()) {
                    std::cout << "error empty" <<endl;
                    return 'f';
                }

                if (!checkSecondaryStructure(ss)) {
                    std::cout << "error bracket" <<endl;
                    return 'n';
                }
            }

            if (!test.compare(keys[1])) {
                std::cout << "sequence: " << it.value() << "\n";
                string seq = it.value();
                if (seq.empty()) {
                    std::cout << "error empty 2" <<endl;
                    return 'l';
                }
            }
            j++;
        }
    std::cout << "no error!\n" << endl;
    }
    return 0;
}*/

bool is_desc_insertible(const string& descfile, const string& rna)
{
    std::ifstream  motif;
    string         line;
    string         seq;
    vector<string> bases;
    int            last;
    char           c    = 'a';
    char*          prev = &c;

    motif = std::ifstream(descfile);
    getline(motif, line);    // ignore "id: number"
    getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    seq  = "";
    last = stoi(bases[1].substr(0, bases[1].find('_')));
    for (vector<string>::iterator b = (bases.begin() + 1); b != (bases.end() - 1); b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = stoi(b->substr(0, b->find('_')));

        if (pos - last > 5) {    // finish this component and start a new one
            seq += ".{5,}";
        } else if (pos - last == 2) {
            seq += ".";
        } else if (pos - last == 3) {
            seq += "..";
        } else if (pos - last == 4) {
            seq += "...";
        } else if (pos - last == 5) {
            seq += "....";
        }
        seq += nt;    // pos - last == 1 in particular
        last = pos;
    }
    smatch m;
    regex  e(seq);

    return regex_search(rna, m, e);

}

vector<vector<Component>> find_next_ones_in(string rna, uint offset, vector<string>& vc, bool is_RIN_or_JSON)
{
    pair<uint, uint>          pos;
    vector<vector<Component>> results;
    vector<vector<Component>> next_ones;
    vector<string>            next_seqs;
    regex                     c(vc[0]);

    //cout << "\t\t>Searching " << vc[0] << " in " << rna << endl;

    if (vc.size() > 1) {
        //cout << "size vc: " << vc.size() << endl; 
        if (regex_search(rna, c)) {
            if (vc.size() > 2) {
                next_seqs = vector<string>(&vc[1], &vc[vc.size()]);
                /*for (uint i = 0; i < next_seqs.size(); i++) {
                    std::cout << "next seq: " << next_seqs[i] << endl;
                }
                std::cout << endl;*/
            }
            else {
                next_seqs = vector<string>(1, vc.back());
                /*for (uint i = 0; i < next_seqs.size(); i++) {
                    std::cout << "next seq: " << next_seqs[i] << endl;
                }
                std::cout << endl;*/
            }
            uint j = 0;
            // For every regexp match
            for (sregex_iterator i = sregex_iterator(rna.begin(), rna.end(), c); i != sregex_iterator(); ++i) {
                smatch match = *i;
                uint add;
                if (is_RIN_or_JSON) add = 2;
                else add = 5;
                pos.first    = match.position() + offset;
                pos.second   = pos.first + match.length() - 1;

                //cout << "\t\t>Inserting " << vc[j] << " in [" << pos.first << ',' << pos.second << "]" << endl;
                // +5 because HL < 3 pbs but not for CaRNAval or Contacts
                // if CaRNAval or Contacts +2 is better
                if (pos.second - offset + add >= rna.length()) {
                     //cout << "\t\t... but we cannot place the next components : Ignored." << endl;
                    continue;
                }
                
                next_ones = find_next_ones_in(rna.substr(pos.second - offset + add), pos.second + add, next_seqs, is_RIN_or_JSON);
                if (!next_ones.size()) {
                    // cout << "\t\t... but we cannot place the next components : Ignored.2" << endl;
                    continue;
                }
                //cout  << endl;
                for (vector<Component> v : next_ones)    // For every combination of the next components
                {
                    // Combine the match for this component pos with the combination
                    // of next_ones as a whole solution
                    vector<Component> r;
                    r.push_back(Component(pos));
                    for (Component& c : v) r.push_back(c);
                    results.push_back(r);
                }
                j++;
            }
        }
    } else {
        // Only one more component to find

        if (regex_search(rna, c)) {
            // For each regexp match
            for (sregex_iterator i = sregex_iterator(rna.begin(), rna.end(), c); i != sregex_iterator(); ++i) {
                smatch match = *i;
                pos.first    = match.position() + offset;
                pos.second   = pos.first + match.length() - 1;

                //cout << "\t\t>Inserting " << vc[0] << " in [" << pos.first << ',' << pos.second << "]" << endl;

                // Create a vector of component with one component for that match
                vector<Component> r;
                r.push_back(Component(pos));
                results.push_back(r);
            }
        }
    }
    return results;
}

bool operator==(const Component& c1, const Component& c2)
{
    if (c1.pos.first != c2.pos.first) return false;
    if (c1.pos.second != c2.pos.second) return false;
    return true;
}

bool operator!=(const Component& c1, const Component& c2) { return not(c1 == c2); }



bool operator==(const Motif& m1, const Motif& m2)
{
    if (m1.get_identifier() != m2.get_identifier()) return false;
    if (m1.score_ != m2.score_) return false;
    if (m1.reversed_ != m2.reversed_) return false;
    for (uint i = 0; i < m1.comp.size(); i++)
        if (m1.comp[i] != m2.comp[i]) return false;
    return true;
}

bool operator!=(const Motif& m1, const Motif& m2) { return not(m1 == m2); }