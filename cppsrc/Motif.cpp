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

    seq  = bases[1].substr(bases[1].find('_') + 1, 1);
    last = std::stoi(bases[1].substr(0, bases[1].find('_')));
    for (vector<string>::iterator b = bases.begin() + 1; b != bases.end() - 1; b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = std::stoi(b->substr(0, b->find('_')));

        if (pos - last > 5) {    // finish this component and start a new one
            seq += nt;
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
        if (pos - last <= 5) seq += nt;
    }
    // Now component_sequences is a vector of sequences like {AGCGC, CGU..GUUU}

    // We need to search for the different positions where to insert the first component
    vector<vector<Component>> vresults = find_next_ones_in(rna, component_sequences);

    // Now create proper motifs
    for (vector<Component>& v : vresults) {
        results.push_back(Motif(v, descfile.substr(0, descfile.find(".desc"))));
    }
    std::cout << "\t>returning vector of size " << results.size() << std::endl;
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

vector<vector<Component>> Motif::find_next_ones_in(string rna, vector<string> vc)
{
    std::smatch               matches;
    std::regex                c(vc[0]);
    pair<uint, uint>          pos;
    vector<vector<Component>> results;
    vector<vector<Component>> next_ones;
    vector<string>            next_seqs;

    if (vc.size() > 1) {
        if (vc.size() > 2)
            next_seqs = vector<string>(&vc[1], &vc[vc.size() - 1]);
        else
            next_seqs = vector<string>(1, vc.back());

        std::regex_search(rna, matches, c);

        for (uint i = 0; i < matches.size(); ++i)    // Pour chacun des matches
        {
            pos.first  = matches.position(i);
            pos.second = matches.length(i) + pos.first - 1;
            std::cout << "\t\t>We can insert " << vc[0] << " in [" << pos.first << ',' << pos.second << ']' << std::endl;
            std::cout << "\t\t>Now searching in " << rna.substr(pos.second + 1) << std::endl;
            next_ones = find_next_ones_in(rna.substr(pos.second + 1), next_seqs);
            for (vector<Component> v : next_ones)    // Pour chacune des combinaisons suivantes
            {
                // Combiner le match et la combinaison suivante
                vector<Component> r;
                r.push_back(Component(pos));
                for (Component& c : v) r.push_back(c);
                results.push_back(r);
            }
        }
    } else {
        std::regex_search(rna, matches, c);

        for (uint i = 0; i < matches.size(); ++i)    // Pour chacun des matches
        {
            pos.first  = matches.position(i);
            pos.second = matches.length(i) + pos.first - 1;
            std::cout << "\t\t>We can insert " << vc[0] << " in [" << pos.first << ',' << pos.second << ']' << std::endl;
            // Combiner le match et la combinaison suivante
            vector<Component> r;
            r.push_back(Component(pos));
            results.push_back(r);
        }
    }
    std::cout << "\t> returning vector of size " << results.size() << std::endl;
    return results;
}

vector<Motif> load_desc_folder(const string& path, const string& rna, bool verbose)
{
    vector<Motif> posInsertionSites;

    if (!exists(path)) {
        std::cerr << "Hmh, i can't find that folder: " << path << std::endl;
        return posInsertionSites;
    } else {
        if (verbose) std::cout << "loading DESC motifs from " << path << "..." << std::endl;
    }

    for (auto it : recursive_directory_range(path)) {
        if (is_desc_insertible(it.path().string(), rna, verbose)) {
            vector<Motif> m = Motif::build_from_desc(it.path().string(), rna);
            for (Motif& mot : m) posInsertionSites.push_back(mot);
        }
    }
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
            std::cout << "\t>Motif " << descfile.substr(0, descfile.find(".desc")) << "   \t" << seq
                      << " can be inserted." << std::endl;
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
