
#ifndef DEF_RNA
#define DEF_RNA

#include <Eigen/Core>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using Eigen::MatrixXf, Eigen::Matrix;
using std::map;
using std::pair;
using std::string;
using std::vector;

#ifndef NUPACK_SHARED_CONSTANTS_H__
enum base_t { BASE_N = 0, BASE_A, BASE_C, BASE_G, BASE_U };    // Comment if you include nupack/shared.h
#else
typedef int base_t;
#endif
enum pair_t { PAIR_AU = 0, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG, PAIR_OTHER = -1 };

class RNA
{
    public:
    RNA(string name, string seq);

    float  get_pij(int i, int j);
    string get_seq(void) const;
    uint   get_RNA_length(void) const;
    void   print_basepair_p_matrix(float theta) const;

    private:
    base_t base_type(char x) const;

    string   name_;    // name of the rna
    string   seq_;     // sequence of the rna with chars
    uint     n_;       // length of the rna
    MatrixXf pij_;     // matrix of basepair probabilities
};

inline float  RNA::get_pij(int i, int j) { return pij_(i, j); }
inline uint   RNA::get_RNA_length() const { return n_; }
inline string RNA::get_seq(void) const { return seq_; }

#endif
