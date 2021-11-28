#ifndef SECONDARY_STRUCTURE_
#define SECONDARY_STRUCTURE_

#include "Motif.h"
#include "rna.h"
#include <iostream>
#include <string>
#include <vector>

using std::pair;
using std::string;
using std::vector;


class SecondaryStructure
{
    public:
    SecondaryStructure(void);
    SecondaryStructure(const RNA& rna);
    SecondaryStructure(bool empty);

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
    size_t        nBP_;           // number of basepairs
    RNA           rna_;           // RNA object which is folded
    bool          is_empty_structure;    // Empty structure, returned when the solver does not find solutions anymore
};

// return if this SecondaryStructure s1 dominates s2
bool operator>(const SecondaryStructure& s1, const SecondaryStructure& s2);
bool operator>=(const SecondaryStructure& s1, const SecondaryStructure& s2);
// return if this SecondaryStructure s2 dominates s1
bool operator<(const SecondaryStructure& s1, const SecondaryStructure& s2);
bool operator<=(const SecondaryStructure& s1, const SecondaryStructure& s2);
// return wether SecondaryStructures are identical or not
bool operator==(const SecondaryStructure& s1, const SecondaryStructure& s2);
bool operator!=(const SecondaryStructure& s1, const SecondaryStructure& s2);

bool motif_sorter(Motif& m1, Motif& m2);
bool basepair_sorter(pair<uint, uint>& i, pair<uint, uint>& j);

inline double SecondaryStructure::get_objective_score(int i) const { return objective_scores_[i - 1]; }
inline void   SecondaryStructure::set_objective_score(int i, double s) { objective_scores_[i - 1] = s; }
inline uint   SecondaryStructure::get_n_motifs(void) const { return motif_info_.size(); }
inline uint   SecondaryStructure::get_n_bp(void) const { return nBP_; }

#endif    //  SECONDARY_STRUCTURE_