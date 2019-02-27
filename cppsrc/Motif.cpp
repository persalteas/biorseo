#include "Motif.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <regex>
#include <sstream>

using namespace boost::filesystem;


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

vector<Motif> Motif::build_from_desc(const string& descfile, string rna)
{
    std::ifstream  motif;
    string         line;
    string         seq;
    vector<string> component_sequences;
    vector<string> bases;
    int            last;
    vector<Motif>  results;
    char           c    = 'a';
    char*          prev = &c;

    motif = std::ifstream(descfile);
    std::getline(motif, line);    // ignore "id: number"
    std::getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    seq  = "";
    last = std::stoi(bases[1].substr(0, bases[1].find('_')));
    for (vector<string>::iterator b = bases.begin() + 1; b != bases.end() - 1; b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = std::stoi(b->substr(0, b->find('_')));

        if (pos - last > 5) {    // finish this component and start a new one
            component_sequences.push_back(seq);
            seq = "";
        } else if (pos - last == 2) {
            seq += '.';
        } else if (pos - last == 3) {
            seq += "..";
        } else if (pos - last == 4) {
            seq += "...";
        } else if (pos - last == 5) {
            seq += "....";
        }
        seq += nt;
        last = pos;
    }
    component_sequences.push_back(seq);
    // Now component_sequences is a vector of sequences like {AGCGC, CGU..GUUU}

    // We need to search for the different positions where to insert the first component
    vector<vector<Component>> vresults = find_next_ones_in(rna, 0, component_sequences);

    // Now create proper motifs
    for (vector<Component>& v : vresults) {
        results.push_back(Motif(v, descfile.substr(0, descfile.find(".desc"))));
    }
    std::cout << results.size() << " times" << std::endl;
    return results;
}

void Motif::load_from_csv(string csv_line)
{
    vector<string> tokens;
    split(tokens, csv_line, boost::is_any_of(","));
    atlas_id = tokens[0];
    score_   = stoi(tokens[2]);
    comp.push_back(Component(std::make_pair<int, int>(stoi(tokens[3]), stoi(tokens[4]))));
    if (tokens[5] != "-") comp.push_back(Component(std::make_pair<int, int>(stoi(tokens[5]), stoi(tokens[6]))));
    reversed_ = (tokens[1] == "True");
    is_model_ = true;
    PDBID     = "";
    source_   = RNAMOTIFATLAS;
}

string Motif::pos_string(void) const
{
    std::stringstream s;
    s << atlas_id << " ( ";
    for (auto c : comp) s << c.pos.first << '-' << c.pos.second << ' ';
    s << ')';
    return s.str();
}

string Motif::get_identifier(void) const
{
    switch (source_) {
    case RNAMOTIFATLAS: return atlas_id; break;
    default: return PDBID;
    }
}

vector<vector<Component>> Motif::find_next_ones_in(string rna, uint offset, vector<string> vc)
{
    pair<uint, uint>          pos;
    vector<vector<Component>> results;
    vector<vector<Component>> next_ones;
    vector<string>            next_seqs;
    std::regex                c(vc[0]);

    // std::cout << "\t\t>Searching " << vc[0] << " in " << rna << std::endl;

    if (vc.size() > 1) {
        if (std::regex_search(rna, c)){
            if (vc.size() > 2)
                next_seqs = vector<string>(&vc[1], &vc[vc.size() - 1]);
            else
                next_seqs = vector<string>(1, vc.back()); 

            // Pour chacun des matches
            for(std::sregex_iterator i = std::sregex_iterator(rna.begin(), rna.end(), c); i != std::sregex_iterator(); ++i )
            {
                std::smatch               match = *i;
                pos.first  = match.position() + offset;
                pos.second = pos.first + match.length() - 1;
                // std::cout << "\t\t>Inserting " << vc[0] << " in [" << pos.first << ',' << pos.second << "]" << std::endl;
                if (pos.second - offset + 5 >= rna.length())
                {
                    // std::cout << "\t\t... but we cannot place the next components : Ignored." << std::endl;
                    continue;
                }
                next_ones = find_next_ones_in(rna.substr(pos.second - offset + 5), pos.second + 5,  next_seqs);
                if (!next_ones.size())
                {
                    // std::cout << "\t\t... but we cannot place the next components : Ignored." << std::endl;
                    continue;
                }
                // std::cout  << std::endl;
                for (vector<Component> v : next_ones)    // Pour chacune des combinaisons suivantes
                {
                    // Combiner le match et la combinaison suivante
                    vector<Component> r;
                    r.push_back(Component(pos));
                    for (Component& c : v) r.push_back(c);
                    results.push_back(r);
                }
            }
        }
    } else {
        if (std::regex_search(rna, c)){
            // Pour chacun des matches
            for(std::sregex_iterator i = std::sregex_iterator(rna.begin(), rna.end(), c); i != std::sregex_iterator(); ++i )
            {
                std::smatch match = *i;
                pos.first  = match.position() + offset;
                pos.second = pos.first + match.length() - 1;
                // std::cout << "\t\t>Inserting " << vc[0] << " in [" << pos.first << ',' << pos.second << "]" << std::endl;
                // Combiner le match et la combinaison suivante
                vector<Component> r;
                r.push_back(Component(pos));
                results.push_back(r);
            }
        }
    }
    return results;
}

char Motif::is_valid_DESC(const string& descfile)
{
    // /!\ returns 0 iff no errors 

    std::ifstream  motif;
    string         line;
    vector<string> bases;
    char           c    = 'a';
    char*          prev = &c;

    motif = std::ifstream(descfile);
    std::getline(motif, line);    // ignore "id: number"
    std::getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    for (vector<string>::iterator b = (bases.begin() + 1); b != (bases.end() - 1); b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = std::stoi(b->substr(0, b->find('_')));

        if (string(1,nt).find_first_not_of("ACGU") != string::npos) return nt;
        if (pos <= 0) return '-';
    }

    while (std::getline(motif, line))
    {
        uint slash = line.find('/');
        string interaction(line.substr(slash-1, 3));

        string b1 = line.substr(line.find('(')+1, 8); // The first base (position_nucleotide)
        string temp = line.substr(slash+1);
        string b2 = temp.substr(temp.find('(')+1, 8); // The second base (position_nucleotide)
        b1.erase(std::remove(b1.begin(), b1.end(), ' '), b1.end());
        b2.erase(std::remove(b2.begin(), b2.end(), ' '), b2.end());
        int p1 = std::stoi(b1.substr(0,b1.find('_')));
        int p2 = std::stoi(b2.substr(0,b2.find('_')));
        
        if ((p2-p1 != 1) and !interaction.compare("C/C")) return 'b';
        if ((p2-p1 <4) and !(interaction.compare("+/+") + interaction.compare("-/-"))) return 'l';
    }
    return 0;
}

vector<Motif> load_desc_folder(const string& path, const string& rna, bool verbose)
{
    vector<Motif> posInsertionSites;
    int errors = 0;
    int accepted = 0;
    int inserted = 0;

    if (!exists(path)) {
        std::cerr << "Hmh, i can't find that folder: " << path << std::endl;
        return posInsertionSites;
    } else {
        if (verbose) std::cout << "loading DESC motifs from " << path << "..." << std::endl;
    }

    char error;
    for (auto it : recursive_directory_range(path)) {
        
        if ((error = Motif::is_valid_DESC(it.path().string()))) {
            std::cerr << "\t>Ignoring motif " << it.path().stem();
            switch (error)
            {
                case '-':
                    std::cerr << ", some nucleotides have a negative number...";
                    break;
                case 'l':
                    std::cerr << ", hairpin (terminal) loops must be at least of size 3 !";
                    break;
                case 'b':
                    std::cerr << ", backbone link between non-consecutive residues ?";
                    break;
                default:
                    std::cerr << ", use of an unknown nucleotide " << error;
            }
            std::cerr << std::endl;
            errors++;
            continue;
        }
        accepted++;
        if (is_desc_insertible(it.path().string(), rna, verbose)) {
            inserted++;
            vector<Motif> m = Motif::build_from_desc(it.path().string(), rna);
            for (Motif& mot : m) posInsertionSites.push_back(mot);
        }
    }
    if (verbose) std::cout << "Inserted " << inserted << " motifs on " << accepted+errors << " (" << errors << " ignored motifs)" << std::endl;
    return posInsertionSites;
}

vector<Motif> load_jar3d_output(const string& path)
{
    vector<Motif> posInsertionSites;
    std::ifstream motifs;
    string        line;

    motifs = std::ifstream(path);
    std::getline(motifs, line);    // skip header
    while (std::getline(motifs, line)) {
        posInsertionSites.push_back(Motif());
        posInsertionSites.back().load_from_csv(line);
    }
    return posInsertionSites;
}

bool is_desc_insertible(const string& descfile, const string& rna, bool verbose)
{
    std::ifstream  motif;
    string         line;
    string         seq;
    vector<string> bases;
    int            last;
    char           c    = 'a';
    char*          prev = &c;

    motif = std::ifstream(descfile);
    std::getline(motif, line);    // ignore "id: number"
    std::getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    seq  = "";
    last = std::stoi(bases[1].substr(0, bases[1].find('_')));
    for (vector<string>::iterator b = (bases.begin() + 1); b != (bases.end() - 1); b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = std::stoi(b->substr(0, b->find('_')));

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
    std::smatch m;
    std::regex  e(seq);
    if (std::regex_search(rna, m, e)) {
        if (verbose)
            std::cout << "\t>Motif " << boost::filesystem::path(descfile).stem() << "   \t" << seq << "\tcan be inserted ";
        return true;
    } else {
        // if (verbose) std::cout << "Ignoring motif " << descfile.substr(0, descfile.find(".desc")) << "   \t" << seq << std::endl;
        return false;
    }
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
