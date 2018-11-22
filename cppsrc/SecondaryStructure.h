#ifndef __INC_IP_SOL__
#define __INC_IP_SOL__

#include "rna.h"
#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector, std::pair;

class SecondaryStructure
{
  public:
    SecondaryStructure(void);
    SecondaryStructure(RNA& rna);

    void   set_basepair(uint i, uint j);
    void   insert_motif(const Motif& m);
    double get_objective_score(int i) const;
    void   set_objective_score(int i, double s);
    void   print(void) const;
    string to_DBN() const;
    string to_string() const;

  private:
    vector<double> objective_scores_;       // values of the different objective functions for that SecondaryStructure
    vector<pair<uint, uint>> basepairs_;    // values of the decision variable of the integer program
    vector<Motif> motif_info_;    // information about known motives in this secondary structure and their positions
    size_t        n_;             // length of the RNA
    size_t        nBP_;
    RNA&          rna_;    // reference to the RNA object which is folded
};

// return if this SecondaryStructure s1 dominates s2
bool operator>(const SecondaryStructure& s1, const SecondaryStructure& s2);
// return if this SecondaryStructure s2 dominates s1
bool operator<(const SecondaryStructure& s1, const SecondaryStructure& s2);

inline double SecondaryStructure::get_objective_score(int i) const { return objective_scores_[i - 1]; }
inline void   SecondaryStructure::set_objective_score(int i, double s) { objective_scores_[i - 1] = s; }

#endif
