#ifndef __INC_IP_SOL__
#define __INC_IP_SOL__

#include "rna.h"
#include <iostream>
#include <string>
#include <vector>

using std::pair;
using std::string;
using std::vector;

typedef struct Comp_ {
    pair<uint, uint> pos;
    int              score;
    size_t           k;
    Comp_(pair<int, int> p, int s) : pos(p), score(s) { k = 1 + pos.second - pos.first; }
} Component;

typedef struct {
    string            atlas_id;
    vector<Component> comp;
    bool              reversed;
    string            pos_string(void) const
    {
        std::stringstream s;
        s << atlas_id << " ( ";
        for (auto c : comp) {
            s << c.pos.first << '-' << c.pos.second << ' ';
        }
        s << ')';
        return s.str();
    }
} Motif;

class SecondaryStructure
{
    public:
    SecondaryStructure(void);
    SecondaryStructure(const RNA& rna);

    void   set_basepair(uint i, uint j);
    void   sort(void);
    void   insert_motif(const Motif& m);
    double get_objective_score(int i) const;
    void   set_objective_score(int i, double s);
    uint   get_n_motifs(void) const;
    uint   get_n_bp(void) const;
    void   print(void) const;
    string to_DBN() const;
    string to_string() const;

    vector<double> objective_scores_;       // values of the different objective functions for that SecondaryStructure
    vector<pair<uint, uint>> basepairs_;    // values of the decision variable of the integer program
    vector<Motif> motif_info_;    // information about known motives in this secondary structure and their positions
    size_t        n_;             // length of the RNA
    size_t        nBP_;
    RNA           rna_;    // RNA object which is folded
};

// return if this SecondaryStructure s1 dominates s2
bool operator>(const SecondaryStructure& s1, const SecondaryStructure& s2);
bool operator>=(const SecondaryStructure& s1, const SecondaryStructure& s2);
// return if this SecondaryStructure s2 dominates s1
bool operator<(const SecondaryStructure& s1, const SecondaryStructure& s2);
bool operator<=(const SecondaryStructure& s1, const SecondaryStructure& s2);
// return wether SecondaryStructures are identical
bool operator==(const SecondaryStructure& s1, const SecondaryStructure& s2);
// utilities to compare secondary structures:
bool operator==(const Motif& m1, const Motif& m2);
bool operator!=(const Motif& m1, const Motif& m2);
bool operator==(const Component& c1, const Component& c2);
bool operator!=(const Component& c1, const Component& c2);
bool motif_sorter(Motif& m1, Motif& m2);
bool basepair_sorter(pair<uint, uint>& i, pair<uint, uint>& j);

inline double SecondaryStructure::get_objective_score(int i) const { return objective_scores_[i - 1]; }
inline void   SecondaryStructure::set_objective_score(int i, double s) { objective_scores_[i - 1] = s; }
inline uint   SecondaryStructure::get_n_motifs(void) const { return motif_info_.size(); }
inline uint   SecondaryStructure::get_n_bp(void) const { return nBP_; }

#endif
