#ifndef MOTIF_H_
#define MOTIF_H_

#include <string>
#include <vector>

using std::pair;
using std::string;
using std::vector;

typedef struct Comp_ {
    pair<uint, uint> pos;
    size_t           k;
    string           seq_;
    Comp_(pair<int, int> p) : pos(p) { k = 1 + pos.second - pos.first; }
} Component;



class Motif
{
    public:
    Motif();
    void              load_from_csv(string csv_line);
    void              build_from_desc(const string& descfile);
    string            pos_string(void) const;
    string            get_origin(void) const;
    string            get_identifier(void) const;
    vector<Component> comp;
    double            score_;
    bool              reversed_;

    private:
    string atlas_id;     // if source = RNAMOTIFATLAS
    string PDBID;        // if source = RNA3DMOTIF
    bool   is_model_;    // Wether the motif is a model or an extracted module from a 3D structure
    enum { RNA3DMOTIF = 1, RNAMOTIFATLAS = 2, CARNAVAL = 3 } source_;
};

bool          is_desc_insertible(const string& descfile, const string& rna);
vector<Motif> load_desc_folder(const string& path);
vector<Motif> load_jar3d_output(const string& path);

// utilities to compare secondary structures:
bool operator==(const Motif& m1, const Motif& m2);
bool operator!=(const Motif& m1, const Motif& m2);
bool operator==(const Component& c1, const Component& c2);
bool operator!=(const Component& c1, const Component& c2);

#endif    // MOTIF_H_