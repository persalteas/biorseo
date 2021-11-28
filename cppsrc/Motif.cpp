#include "Motif.h"
#include "Pool.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <regex>
#include <sstream>
#include <thread>

using namespace boost::filesystem;
using namespace std;
using json = nlohmann::json;

uint Motif::delay = 1;

Motif::Motif(void) {}

Motif::Motif(string csv_line)
{
    vector<string> tokens;
    split(tokens, csv_line, boost::is_any_of(","));

    id_ = tokens[0];
    score_  = stoi(tokens[1]);
    source_ = CSV;
    reversed_ = false; // default value

    if (csv_line.find(string("True")) != std::string::npos or csv_line.find(string("False")) != std::string::npos) 
    {
        // This has been created by jar3d
        score_   = stoi(tokens[2]);
        comp.push_back(Component(make_pair<int, int>(stoi(tokens[3]), stoi(tokens[4]))));
        if (tokens[5] != "-") comp.push_back(Component(make_pair<int, int>(stoi(tokens[5]), stoi(tokens[6]))));
        reversed_ = (tokens[1] == "True");
    }
    else // this has been created by BayesPairing
    {
        uint i = 2;
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

Motif::Motif(const vector<Component>& v, string name) : comp(v), id_(name)
{
    source_ = RNA3DMOTIF;
    reversed_ = false;
}

Motif::Motif(const vector<Component>& v, string name, string& struc) : comp(v), id_(name)
{
    source_ = JSON;
    reversed_ = false;
    
    stack<uint> rbrackets;
    stack<uint> sbrackets;
    stack<uint> cbrackets;
    stack<uint> chevrons;
    
    uint count = 0;
    uint debut = v[count].pos.first;
    uint gap = 0;

    for (uint i = 0; i < struc.size(); i++) {
        if (struc[i] == '(') {
            rbrackets.push(i + debut + gap - count);
        } else if (struc[i] == ')') {
            Link l;
            l.nts.first = rbrackets.top();
            l.nts.second = i + debut + gap - count;
            links_.push_back(l);
            rbrackets.pop();
        } else if (struc[i] == '[') {
            sbrackets.push(i + debut + gap - count);
        } else if (struc[i] == ']') {
            Link l;
            l.nts.first = sbrackets.top();
            l.nts.second = i + debut + gap - count;
            links_.push_back(l);
            sbrackets.pop();
        } else if (struc[i] == '{') {
            cbrackets.push(i + debut + gap - count);
        } else if (struc[i] == '}') {
            Link l;
            l.nts.first = cbrackets.top();
            l.nts.second = i + debut + gap - count;
            links_.push_back(l);
            cbrackets.pop();
        } else if (struc[i] == '<') {
            chevrons.push(i + debut + gap - count);
        } else if (struc[i] == '>') {
            Link l;
            l.nts.first = chevrons.top();
            l.nts.second = i + debut + gap - count;
            links_.push_back(l);
            chevrons.pop();
        } else if (struc[i] == '&') {
            count ++;
            gap += v[count].pos.first - v[count - 1].pos.second - 1;
       }
    }
}

Motif::Motif(const vector<Component>& v, path rinfile, uint id, bool reversed) : comp(v), reversed_(reversed)
{
    // Loads a motif from the RIN file of Carnaval
    id_     = to_string(id);
    source_ = CARNAVAL;

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

        // Now renumber the Links based on the components positions.
        for (Link& l : links_) {
            size_t sum_comp = 0;
            size_t j;
            for (j=0; j<v.size(); j++) {
                const Component& c = v[j];
                if (l.nts.first >= sum_comp and l.nts.first < sum_comp + c.k) {
                    // This is the right component
                    l.nts.first += c.pos.first - sum_comp;
                    sum_comp += c.k;
                    break;
                }
                sum_comp += c.k;
            }

            for (; j<v.size(); j++) {
                const Component& c = v[j];
                if (l.nts.second >= sum_comp and l.nts.second < sum_comp + c.k) {
                    // This is the right component
                    l.nts.second += c.pos.second - sum_comp;
                    break;
                }
                sum_comp += c.k;
            }
        }
    }
    else cout << "\t> RIN file not found : " << rinfile << endl;
}

string Motif::pos_string(void) const
{
    stringstream s;
    s << id_ << " ( ";
    for (auto c : comp) s << c.pos.first << '-' << c.pos.second << ' ';
    s << ')';
    return s.str();
}

string Motif::sec_struct(void) const 
{
    // Fill a string of the right size with dots
    string secstruct = string("");
    secstruct += to_string(comp.size()) + string(" components: ");
    for (auto c : comp) {
        for (uint k = c.pos.first; k<=c.pos.second; k++) secstruct += ".";
        secstruct += " ";
    }

    // Replace dots by brackets
    secstruct += "basepairs:";
    for (auto l : links_) {
        secstruct += '\t' + to_string(l.nts.first) + '-' + to_string(l.nts.second);
    }
    return secstruct;
}

string Motif::get_identifier(void) const
{
    switch (source_) {
    case CARNAVAL: return string("RIN") + id_; break;
    case JSON: return string("JSON") + id_; break;
    default: return id_;
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
        // components are formatted like:
        // pos;k;seq (position, length, sequence)
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

char Motif::is_valid_JSON(const json_elem& i)
{
    string seq = i.value()["sequence"];
    string ss = i.value()["struct2d"];

    if (ss.size() != seq.size()) return 'x';
    if (!check_motif_sequence(seq)) return 'r';
    if (!check_motif_ss(ss)) return 'n';
    if (seq.empty()) return 'e';
    if (seq.size() - count(seq.begin(),seq.end(), '&') < 4) return 'l';
    
    // Iterate on components to check their length
    string subseq;
    vector<string> components;
    uint end, n=0;
    while((seq.find('&') != string::npos)) {
        end = seq.find('&');  
        subseq = seq.substr(0, end);
        seq = seq.substr(end + 1);
        if (subseq.size() >= 2) components.push_back(subseq); 
    }
    if (seq.size() >= 2) // Last component after the last &
        components.push_back(seq); 
    for (auto comp : components)
        n += comp.size();
    if(n <= 3) return 'l';

    return char(0);
}

//Check if the sequence is a rna sequence (ATUGC-only)
bool check_motif_sequence(string seq) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    for (char c : seq)
        if (!(c == 'A' || c == 'U' || c == 'T' || c == '&' || c == 'G' || c == 'C')) {
            return false;
        }
    return true;
}

// Check that there are as many opening brackets as closing ones
bool check_motif_ss(string struc)
{ 
    stack<uint> rbrackets;
    stack<uint> sbrackets;
    stack<uint> cbrackets;
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
                    rbrackets.push(i);

                } else if (struc[i] == ')') {
                    if (!rbrackets.empty())
                        rbrackets.pop();
                    else return false;

                } else if (struc[i] == '[') {
                    sbrackets.push(i);

                } else if (struc[i] == ']') {
                    if (!sbrackets.empty())
                        sbrackets.pop();
                    else return false;

                } else if (struc[i] == '{') {
                    cbrackets.push(i);

                } else if (struc[i] == '}') {
                    if (!cbrackets.empty())
                        cbrackets.pop();
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
    return (rbrackets.empty() && sbrackets.empty() && cbrackets.empty() && chevrons.empty());
}

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

vector<vector<Component>> find_next_ones_in(string rna, uint offset, vector<string> vc)
{
    pair<uint, uint>          pos;
    vector<vector<Component>> results;
    vector<vector<Component>> next_ones;
    vector<string>            next_seqs;
    regex                     c(vc[0]);

    if (vc.size() > 1) {
        if (regex_search(rna, c)) {

            // Define a vector with the remaining components to insert behind vc[0]
            if (vc.size() > 2) {
                next_seqs = vector<string>(&vc[1], &vc[vc.size()]);
            }
            else {
                next_seqs = vector<string>(1, vc.back());
            }

            // For every regexp match
            for (sregex_iterator i = sregex_iterator(rna.begin(), rna.end(), c); i != sregex_iterator(); ++i) {
                //retrieve the hit position
                smatch match = *i;
                pos.first    = match.position() + offset;
                pos.second   = pos.first + match.length() - 1;

                // Check if we can insert the next components after this hit for vc[0]
                // +5 because HL < 3 pbs but not for CaRNAval or Contacts
                // if CaRNAval or Contacts +2 is better
                if (pos.second - offset + Motif::delay >= rna.length()) break; // no room for the next components. give up.
                next_ones = find_next_ones_in(rna.substr(pos.second - offset + Motif::delay), pos.second + Motif::delay, next_seqs);
                if (!next_ones.size()) break; // no matches for the next components. give up.

                 // For every solution for the next components, combine this solution for vc[0] (pos) with the one for the following vc[i] (v)
                for (vector<Component> v : next_ones)
                {
                    // Combine the match for this component pos with the combination
                    // of next_ones as a whole solution
                    vector<Component> r;
                    r.push_back(Component(pos));    // solution for vc[0]
                    for (Component& c : v) r.push_back(c); // solutions for the following vc[i]
                    results.push_back(r);
                }
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