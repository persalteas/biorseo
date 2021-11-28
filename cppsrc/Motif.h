#ifndef MOTIF_H_
#define MOTIF_H_

#include <boost/filesystem.hpp>
#include <mutex>
#include <string>
#include <vector>
#include <filesystem>
#include "rna.h"
#include "json.hpp"

using boost::filesystem::path;
using std::pair;
using std::string;
using std::vector;
using std::mutex;

typedef enum { RNA3DMOTIF = 1, CSV = 2, CARNAVAL = 3, JSON = 4 } source_type;
typedef nlohmann::detail::iter_impl<nlohmann::basic_json<> > json_elem;


typedef struct Comp_ {
    pair<uint, uint> pos;
    size_t           k;
    string           seq_;
    Comp_(pair<int, int> p) : pos(p) { k = 1 + pos.second - pos.first; }
    Comp_(uint start, uint length) : k(length)
    {
        pos.first  = start;
        pos.second = start + length - 1;
    }
} Component;


typedef struct Link
{
    pair<uint, uint> nts;
    bool long_range;
} Link ;



class Motif
{
  public:
    Motif(void);
    Motif(string csv_line);
    Motif(const vector<Component>& v, string name);
    Motif(const vector<Component>& v, string name, string& struc);
    Motif(const vector<Component>& v, path rinfile, uint id, bool reversed);
    // Motif(string path, int id); //full path to biorseo/data/modules/RIN/Subfiles/
    
    static char       is_valid_RIN(const string& rinfile);
    static char       is_valid_DESC(const string& descfile);
    static char       is_valid_JSON(const json_elem& i);

    string            pos_string(void) const;
    string            sec_struct(void) const;
    string            get_origin(void) const;
    string            get_identifier(void) const;
    vector<Component> comp;
    vector<Link>      links_;
    vector<uint>      pos_contacts;

    size_t            contact_;
    double            tx_occurrences_;
    double            score_;
    bool              reversed_;
    static uint       delay;
    // delay is the minimal shift between end of a component and begining of the next.
    // For regular loop motifs, it should be at least 5 (because hairpins cannot be of size smaller than 5).
    // For the general case, it could be zero, but solutions will look dirty...
    // Higher values reduce combinatorial explosion of potential insertion sites.

  private:
    string id_;
    source_type source_;
};

bool                        is_desc_insertible(const string& descfile, const string& rna);
bool                        check_motif_ss(string);
bool                        check_motif_sequence(string);

vector<Motif>               load_txt_folder(const string& path, const string& rna, bool verbose);
vector<Motif>               load_desc_folder(const string& path, const string& rna, bool verbose);
vector<Motif>               load_csv(const string& path);
vector<Motif>               load_json_folder(const string& path, const string& rna, bool verbose);

vector<vector<Component>>   find_next_ones_in(string rna, uint offset, vector<string> vc);

// utilities to compare secondary structures:
bool operator==(const Motif& m1, const Motif& m2);
bool operator!=(const Motif& m1, const Motif& m2);
bool operator==(const Component& c1, const Component& c2);
bool operator!=(const Component& c1, const Component& c2);

#endif    // MOTIF_H_