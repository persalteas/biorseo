#include "Motif.h"
#include "Pool.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <regex>
#include <sstream>
#include <thread>

using namespace boost::filesystem;
using namespace std;


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

void Motif::build_from_desc(args_of_parallel_func arg_struct)
// void Motif::build_from_desc(path descfile, string rna, vector<Motif>& final_results)
{
    path           descfile                 = arg_struct.descfile;
    string&        rna                      = arg_struct.rna;
    vector<Motif>& final_results            = arg_struct.final_results;
    mutex&         posInsertionSites_access = arg_struct.posInsertionSites_mutex;

    std::ifstream             motif;
    vector<vector<Component>> vresults;
    string                    line;
    string                    seq;
    vector<string>            component_sequences;
    vector<string>            bases;
    int                       last;
    char                      c    = 'a';
    char*                     prev = &c;

    motif = std::ifstream(descfile.string());
    getline(motif, line);    // ignore "id: number"
    getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    seq  = "";
    last = stoi(bases[1].substr(0, bases[1].find('_')));
    for (vector<string>::iterator b = bases.begin() + 1; b != bases.end() - 1; b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = stoi(b->substr(0, b->find('_')));

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

    // identify components of length 1 or 2 to extend them to length 3
    vector<uint> comp_of_size_1;
    vector<uint> comp_of_size_2;
    for (uint p = 0; p < component_sequences.size(); ++p) {
        if (component_sequences[p].length() == 1) comp_of_size_1.push_back(p);
        if (component_sequences[p].length() == 2) comp_of_size_2.push_back(p);
    }
    if (comp_of_size_1.size() or comp_of_size_2.size()) {
        vector<vector<string>> motif_variants;    // Will contain several component_sequences vectors according to the size where you extend too short components

        component_sequences.clear();    // rebuild from scratch
        motif_variants.push_back(component_sequences);
        uint actual_comp = 0;

        seq  = "";
        last = stoi(bases[1].substr(0, bases[1].find('_')));
        for (vector<string>::iterator b = bases.begin() + 1; b < bases.end() - 1; b++) {
            int  pos = stoi(b->substr(0, b->find('_')));
            char nt  = b->substr(b->find('_') + 1, 1).back();
            if (comp_of_size_1.size() and actual_comp == comp_of_size_1[0])    // we are on the first component of size 1
            {
                b--;
                nt          = b->substr(b->find('_') + 1, 1).back();
                string seq1 = "";
                seq1 += nt;
                seq1 += "..";
                string seq2 = ".";
                seq2 += nt;
                seq2 += ".";
                string seq3 = "..";
                seq3 += nt;
                uint end = motif_variants.size();    // before to add the new ones
                for (uint u = 0; u < end; ++u) {
                    motif_variants.push_back(motif_variants[u]);    // copy 1 for seq2
                    motif_variants.back().push_back(seq2);
                    motif_variants.push_back(motif_variants[u]);    // copy 2 for seq3
                    motif_variants.back().push_back(seq3);
                    motif_variants[u].push_back(seq1);
                }
                seq = "";
                actual_comp++;
                comp_of_size_1.erase(comp_of_size_1.begin());    // the first element has been processed, remove it
                last = pos;
            } else if (comp_of_size_2.size() and actual_comp == comp_of_size_2[0]) {    // we are on the first component of size 2
                b--;
                nt = b->substr(b->find('_') + 1, 1).back();
                b++;    // skip the next nucleotide
                char next   = b->substr(b->find('_') + 1, 1).back();
                last        = stoi(b->substr(0, b->find('_')));
                string seq1 = "";
                seq1 += nt;
                seq1 += next;
                seq1 += ".";
                string seq2 = ".";
                seq2 += nt;
                seq2 += next;
                uint end = motif_variants.size();    // before to add the new one
                for (uint u = 0; u < end; ++u) {
                    motif_variants.push_back(motif_variants[u]);    // copy 1 for seq2
                    motif_variants.back().push_back(seq2);
                    motif_variants[u].push_back(seq1);
                }
                seq = "";
                actual_comp++;
                comp_of_size_2.erase(comp_of_size_2.begin());    // the first element has been processed, remove it
            } else {                                             // we are on a longer component
                if (pos - last > 5) {                            // finish this component and start a new one
                    actual_comp++;
                    for (vector<string>& c_s : motif_variants) c_s.push_back(seq);
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
        }
        for (auto c_s : motif_variants)
            if (seq.length()) c_s.push_back(seq);    // pushing the last one after iterating over the bases

        // We need to search for the different positions where to insert the first component
        for (auto c_s : motif_variants) {
            vector<vector<Component>> new_results = find_next_ones_in(rna, 0, c_s);
            vresults.insert(vresults.end(), new_results.begin(), new_results.end());
        }

    } else {
        // No multiple motif variants : we serach in a single vector component_sequences
        // We need to search for the different positions where to insert the first component
        vresults = find_next_ones_in(rna, 0, component_sequences);
    }

    // Now create proper motifs with Motif class
    for (vector<Component>& v : vresults) {
        unique_lock<mutex> lock(posInsertionSites_access);
        final_results.push_back(Motif(v, path(descfile).stem().string()));
        lock.unlock();
    }
}

void Motif::load_from_csv(string csv_line)
{
    vector<string> tokens;
    split(tokens, csv_line, boost::is_any_of(","));
    if (csv_line.find(string("True")) != std::string::npos or csv_line.find(string("False")) != std::string::npos) {    // This has been created by jar3d
        atlas_id = tokens[0];
        score_   = stoi(tokens[2]);
        comp.push_back(Component(make_pair<int, int>(stoi(tokens[3]), stoi(tokens[4]))));
        if (tokens[5] != "-") comp.push_back(Component(make_pair<int, int>(stoi(tokens[5]), stoi(tokens[6]))));
        reversed_ = (tokens[1] == "True");
        is_model_ = true;
        PDBID     = "";
        source_   = RNAMOTIFATLAS;
    } else {    // this has been created by BayesPairing
        score_ = stoi(tokens[1]);
        // identify source:
        if (tokens[0].find(string("rna3dmotif")) == std::string::npos) {
            is_model_ = true;
            PDBID     = "";
            source_   = RNAMOTIFATLAS;
            atlas_id  = tokens[0];
        } else {
            is_model_ = false;
            PDBID     = tokens[0];
            source_   = RNA3DMOTIF;
            atlas_id  = "";
        }
        uint i = 2;
        while (i < tokens.size()) {
            if (stoi(tokens[i]) < stoi(tokens[i + 1]))
                comp.push_back(Component(make_pair<int, int>(stoi(tokens[i]), stoi(tokens[i + 1]))));
            else
                i = i + 0;
            i += 2;
        }
    }
}



void Motif::load_from_txt(string path, int id)
{
    carnaval_id = to_string(id) ;
    atlas_id = "" ;
    PDBID = "" ;
    is_model_ = true ;


    /*-----comp-----*/
    std::ifstream file(path + carnaval_id + ".txt") ;

    if (file.is_open())
    {
        string line ;
        std::getline(file,line) ; //skip the header line

        size_t index = 0 ;

        string pos_str ;
            size_t sub_index = 0 ;
            string sub_pos_str ;

        string k_str ;

        string seq ;

        while ( std::getline(file,line) )
        {
            if (line == "\n") break ; //skip last line (empty)

            Component c(0,0) ;

            //c.pos
            index = line.find(";") ;
            pos_str = line.substr(0, index) ;
            line.erase(0, index+1) ;

                sub_index = pos_str.find(",") ;
                sub_pos_str = pos_str.substr(0, sub_index) ;
                pos_str.erase(0, sub_index+1) ;
                c.pos.first = std::stoi(sub_pos_str) ;
                c.pos.second = std::stoi(pos_str) ;

            //c.k
            index = line.find(";") ;
            k_str = line.substr(0, index) ;
            line.erase(0, index+1) ;
            c.k = std::stoi(k_str) ;

            //c.seq_
            seq = line ;
            c.seq_ = seq ;

            comp.push_back(c) ;
        }
    }

    else std::cout << "Motif::load_from_txt -> File not found : " + path + carnaval_id + ".txt" << std::endl ;
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
    case CARNAVAL return carnaval_id; break;
    default: return PDBID;
    }
}

vector<vector<Component>> Motif::find_next_ones_in(string rna, uint offset, vector<string> vc)
{
    pair<uint, uint>          pos;
    vector<vector<Component>> results;
    vector<vector<Component>> next_ones;
    vector<string>            next_seqs;
    regex                     c(vc[0]);

    // cout << "\t\t>Searching " << vc[0] << " in " << rna << endl;

    if (vc.size() > 1) {
        if (regex_search(rna, c)) {
            if (vc.size() > 2)
                next_seqs = vector<string>(&vc[1], &vc[vc.size() - 1]);
            else
                next_seqs = vector<string>(1, vc.back());

            // Pour chacun des matches
            for (sregex_iterator i = sregex_iterator(rna.begin(), rna.end(), c); i != sregex_iterator(); ++i) {
                smatch match = *i;
                pos.first    = match.position() + offset;
                pos.second   = pos.first + match.length() - 1;
                // cout << "\t\t>Inserting " << vc[0] << " in [" << pos.first << ',' << pos.second << "]" << endl;
                if (pos.second - offset + 5 >= rna.length()) {
                    // cout << "\t\t... but we cannot place the next components : Ignored." << endl;
                    continue;
                }
                next_ones = find_next_ones_in(rna.substr(pos.second - offset + 5), pos.second + 5, next_seqs);
                if (!next_ones.size()) {
                    // cout << "\t\t... but we cannot place the next components : Ignored." << endl;
                    continue;
                }
                // cout  << endl;
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
        if (regex_search(rna, c)) {
            // Pour chacun des matches
            for (sregex_iterator i = sregex_iterator(rna.begin(), rna.end(), c); i != sregex_iterator(); ++i) {
                smatch match = *i;
                pos.first    = match.position() + offset;
                pos.second   = pos.first + match.length() - 1;
                // cout << "\t\t>Inserting " << vc[0] << " in [" << pos.first << ',' << pos.second << "]" << endl;
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



vector<Motif> Motif::RIN_list(const string& rna, bool reversed)
{
    string used_rna = rna ;
    if (reversed) std::reverse(used_rna.begin(), used_rna.end()) ;

    vector<Motif> res;
    vector< vector<int> > comps_starts ; //list of beginning indexes by components
    vector<int> ks;
    size_t index ;

    for (Component& component : comp)
    {
        index = 0 ;
        vector<int> starts ;

        while (index != string::npos)
        {
            index = used_rna.find(component.seq_, index) ;
            starts.push_back(index);
        }

        comps_starts.push_back(starts) ;
        ks.push_back(component.k) ;
    }


    return res ;
}



bool Motif::is_valid(const string& rna, bool reversed) //renvoyer un vecteur de motifs
{
    size_t index = 0 ;
    reversed_ = reversed ;

    for (Component& component : comp)
    {
        index = rna.find(component.seq_, index) ;

        if (index == string::npos) //seq_ not found
        {
            if (reversed) return false ; //searched through rna and reversed rna, still nothing
            else //try to look through the reversed sequence
            {
                string reversed_rna(rna) ;
                std::reverse(reversed_rna.begin(), reversed_rna.end()) ;
                return is_valid(reversed_rna, true) ;
            }
        }

        component.pos.first = index ;
        component.pos.second = index + component.k - 1 ;
    }

    return true ;
}



vector<Motif> load_txt_folder(const string& path, const string& rna, bool verbose)
{
    vector<Motif> motifs;
    string valid_path = path ;

    string reversed_rna = rna ;
    std::reverse(reversed_rna.begin(), reversed_rna.end()) ;

    bool verified ;
    

    if (valid_path.back() != '/') valid_path.push_back('/') ;

    if (!exists(valid_path))
    {
        cerr << "Hmh, i can't find that folder: " << valid_path << endl;
        return motifs;
    }
    else if (verbose) cout << "loading RIN motifs from " << valid_path << "..." << endl;

    for (int i=0; i<337; i++) //337 is the number of RINs in CaRNAval
    {
        motifs.push_back(Motif()) ;
        motifs.back().load_from_txt(valid_path, i);

        vector<string> vc; 

        for (Component component : motifs.back().comp)
            vc.push_back(component.seq_) ;

        vector<vector<Component>> occurrences = motifs.back().find_next_ones_in(rna, 0, vc) ;
        vector<vector<Component>> r_occurrences = motifs.back().find_next_ones_in(reversed_rna, 0, vc) ;
        motifs.pop_back() ;

        for (vector<Component> occ : occurrences)
        {
            motifs.push_back(Motif()) ;
            motifs.back().load_from_txt(valid_path, i);
            motifs.back().comp = occ ;
            motifs.back().reversed_ = false ;
        }

        for (vector<Component> occ : r_occurrences)
        {
            motifs.push_back(Motif()) ;
            motifs.back().load_from_txt(valid_path, i);
            motifs.back().comp = occ ;
            motifs.back().reversed_ = true ;
        }


        /*
        verified = motifs.back().is_valid(rna, false) ;
        if ( !verified ) //if the motif is not in the RNA sequence, we remove it
        {
            motifs.pop_back() ;
            if (verbose) std::cout << "RIN n°" << i+1 << " not found in RNA sequence\n" ;
        }
        else if (verbose)
        {
            std::cout << "RIN n°" << i+1 << " found at : " ;
            for (Component component : motifs.back().comp)
            {
                std::cout << component.pos.first << "-" << component.pos.second << " " ;
            }
            std::cout << std::endl ;
        }*/
    }

    if (verbose) cout << "Done" << endl;
    
    return motifs ;
}



vector<Motif> load_desc_folder(const string& path, const string& rna, bool verbose)
{
    vector<Motif> posInsertionSites;
    mutex         posInsertionSites_access;
    Pool          pool;
    int           errors   = 0;
    int           accepted = 0;
    int           inserted = 0;
    // int           num_threads = thread::hardware_concurrency();
    int            num_threads = 2;
    vector<thread> thread_pool;

    if (!exists(path)) {
        cerr << "Hmh, i can't find that folder: " << path << endl;
        return posInsertionSites;
    } else {
        if (verbose) cout << "loading DESC motifs from " << path << "..." << endl;
    }

    for (int i = 0; i < num_threads; i++) thread_pool.push_back(thread(&Pool::infinite_loop_func, &pool));

    char error;
    for (auto it : recursive_directory_range(path)) {    // Add every .desc file to the queue (iff valid)
        if ((error = Motif::is_valid_DESC(it.path().string()))) {
            if (verbose) {
                cerr << "\t>Ignoring motif " << it.path().stem();
                switch (error) {
                case '-': cerr << ", some nucleotides have a negative number..."; break;
                case 'l': cerr << ", hairpin (terminal) loops must be at least of size 3 !"; break;
                case 'b': cerr << ", backbone link between non-consecutive residues ?"; break;
                default: cerr << ", use of an unknown nucleotide " << error;
                }
                cerr << endl;
            }
            errors++;
            continue;
        }
        accepted++;
        if (is_desc_insertible(it.path().string(), rna, verbose)) {
            args_of_parallel_func args(it.path(), rna, posInsertionSites, posInsertionSites_access);
            inserted++;
            pool.push(bind(Motif::build_from_desc, args));
            // Motif::build_from_desc(it.path(), rna, posInsertionSites);
        }
    }
    pool.done();
    for (unsigned int i = 0; i < thread_pool.size(); i++) thread_pool.at(i).join();
    if (verbose)
        cout << "Inserted " << inserted << " motifs on " << accepted + errors << " (" << errors << " ignored motifs)" << endl;
    return posInsertionSites;
}

vector<Motif> load_csv(const string& path)
{
    vector<Motif> posInsertionSites;
    std::ifstream motifs;
    string        line;

    motifs = std::ifstream(path);
    getline(motifs, line);    // skip header
    while (getline(motifs, line)) {
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
    if (regex_search(rna, m, e)) {
        if (verbose)
            cout << "\t>Motif " << boost::filesystem::path(descfile).stem() << "   \t" << seq << "\tcan be inserted " << endl;
        return true;
    } else {
        // if (verbose) cout << "Ignoring motif " << descfile.substr(0, descfile.find(".desc")) << "   \t" << seq << endl;
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
