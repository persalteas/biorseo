#ifndef MOTIF_H_
#define MOTIF_H_

#include <boost/filesystem.hpp>
#include <mutex>
#include <string>
#include <vector>

using boost::filesystem::path;
using std::pair;
using std::string;
using std::vector;

class Motif;    // forward declaration

typedef struct args_ {
    path           descfile;
    string         rna;
    vector<Motif>& final_results;
    std::mutex&    posInsertionSites_mutex;
    args_(path descfile_, string rna_, vector<Motif>& vector_, std::mutex& mutex_)
    : descfile(descfile_), rna(rna_), final_results(vector_), posInsertionSites_mutex(mutex_)
    {
    }
} args_of_parallel_func;

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

class Motif
{
    public:
    Motif();
    Motif(const vector<Component>& v, string PDB);
    void load_from_json(int id);
    void load_from_csv(string csv_line);
    // static void       build_from_desc(path descfile, string rna, vector<Motif>& final_results);
    static void       build_from_desc(args_of_parallel_func args);
    static char       is_valid_DESC(const string& descfile);
    string            pos_string(void) const;
    string            get_origin(void) const;
    string            get_identifier(void) const;
    vector<Component> comp;
    double            score_;
    bool              reversed_;

    private:
    static vector<vector<Component>> find_next_ones_in(string rna, uint offset, vector<string> vc);

    string carnaval_id;  // if source = CARNAVAL
    string atlas_id;     // if source = RNAMOTIFATLAS
    string PDBID;        // if source = RNA3DMOTIF
    bool   is_model_;    // Wether the motif is a model or an extracted module from a 3D structure
    enum { RNA3DMOTIF = 1, RNAMOTIFATLAS = 2, CARNAVAL = 3 } source_;
};

bool          is_desc_insertible(const string& descfile, const string& rna, bool verbose);
vector<Motif> load_desc_folder(const string& path, const string& rna, bool verbose);
vector<Motif> load_csv(const string& path);

// utilities to compare secondary structures:
bool operator==(const Motif& m1, const Motif& m2);
bool operator!=(const Motif& m1, const Motif& m2);
bool operator==(const Component& c1, const Component& c2);
bool operator!=(const Component& c1, const Component& c2);

#endif    // MOTIF_H_