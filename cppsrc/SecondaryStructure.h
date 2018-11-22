#ifndef __INC_IP_SOL__
#define __INC_IP_SOL__

#include "rna.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

typedef vector<int>                 VI;
typedef vector<VI>                  VVI;
typedef vector<std::pair<int, int>> VII;


class SecondaryStructure
{
  public:
    SecondaryStructure(void);
    SecondaryStructure(const vector<double>& scores, const vector<bool>& decision_variables, VII coord, int RNAlength);

    double              get_objective_score(int i) const;
    const vector<bool>& get_decision_variables() const;
    bool                get_decision_value(int i) const;
    VII                 get_coord() const;
    int                 get_RNA_length() const;

    void   set_objective_score(int i, double s);
    string to_DBN() const;
    string to_string() const;

  private:
    vector<double> objective_scores_;    // values of the different objective functions for that SecondaryStructure
    vector<bool>   dv_;                  // values of the decision variable of the integer program
    vector<Motif>  motif_info_;    // information about known motives in this secondary structure and their positions
    VII            coord_;         // coordinates of the dv_. dv_[i] == true <==> coord_[i][0] paired to coord_[i][1];
    size_t         n_;             // length of the RNA
};

// return if this SecondaryStructure s1 dominates s2
bool operator>(const SecondaryStructure& s1, const SecondaryStructure& s2);
// return if this SecondaryStructure s2 dominates s1
bool operator<(const SecondaryStructure& s1, const SecondaryStructure& s2);

inline double              SecondaryStructure::get_objective_score(int i) const { return objective_scores_[i]; }
inline const vector<bool>& SecondaryStructure::get_decision_variables() const { return dv_; }
inline void                SecondaryStructure::set_objective_score(int i, double s) { objective_scores_[i - 1] = s; }
inline VII                 SecondaryStructure::get_coord() const { return coord_; }
inline int                 SecondaryStructure::get_RNA_length() const { return n_; }
inline bool                SecondaryStructure::get_decision_value(int i) const { return dv_[i]; }

#endif
