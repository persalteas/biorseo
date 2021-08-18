
#ifndef DEF_RNA
#define DEF_RNA

#include <eigen3/Eigen/Core>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using Eigen::MatrixXf;
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
    RNA(void);
    RNA(string name, string seq, bool verbose);

    float  get_pij(int i, int j);
    string get_seq(void) const;
    uint   get_RNA_length(void) const;
    void   print_basepair_p_matrix(float theta) const;

    vector<pair<int,int>> get_coord();
    vector<vector<int>> get_type();
    int find_coord(pair<int, int>);

    bool verbose_;    // Should we print things ?

    private:
    base_t base_type(char x) const;

    string   name_;    // name of the rna
    string   seq_;     // sequence of the rna with chars
    uint     n_;       // length of the rna
    MatrixXf pij_;     // matrix of basepair probabilities

    vector<vector<int>> type_;  //vector of base pair types
    vector<pair<int,int>> coord_; //vector of base pair coordinates
};

inline float  RNA::get_pij(int i, int j) { return pij_(i, j); }
inline uint   RNA::get_RNA_length() const { return n_; }
inline string RNA::get_seq(void) const { return seq_; }

inline vector<pair<int,int>>  RNA::get_coord() { return coord_; }
inline vector<vector<int>>  RNA::get_type() { return type_; }
inline int RNA::find_coord(std::pair< int,int > p){
    std::vector<pair<int, int>>::iterator it = find(coord_.begin(), coord_.end(), p);
    int r = -1;
    if(it != coord_.end())
        r = std::distance(coord_.begin(), it);
    return r;
}

#endif
