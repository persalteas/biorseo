
#ifndef DEF_RNA
#define DEF_RNA

#include <Eigen/CXX11/Tensor>
#include <Eigen/Core>
#include <boost/multi_array.hpp>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#define kB 1.380e-23               // Boltzmann constant in kcal/mol/K
#define ZERO_C_IN_KELVIN 273.15    // Zero degrees C in Kelvin
#define AVOGADRO 6.022e23          // Avogadro's number

using boost::multi_array;
using Eigen::MatrixXf;

using std::map;
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
        for (auto c : comp) {
            s << c.pos.first << '-' << c.pos.second << ' ';
        }
        return s.str();
    }
} Motif;

typedef struct {
    float a1;    // multiloop penalty
    float a2;    // penalty of closing basepairs in the multiploop
    float a3;    // penalty for unpaired nt in the multiloop
    float hairpin37[30];
    float bulge37[30];
    float interior37[30];
    float stack37[6][6];
    float int11_37[6][6][4][4];
    float int21_37[6][4][4][6][4];
    float int22_37[6][6][4][4][4][4];
    float dangle3_37[6][4];
    float dangle5_37[6][4];
    float triloop37[4][4][4][4][4];
    float tloop37[4][4][4][4][4][4];
    float mismatch_hairpin37[4][4][6];
    float mismatch_interior37[4][4][6];
    float asymmetry_penalty[4];
    float polyC_penalty, polyC_slope, polyC_int;
    float at_penalty;
    float pk_penalty;              // beta1
    float pk_multiloop_penalty;    // beta1m
    float pk_pk_penalty;           // beta1p
    float pk_paired_penalty;       // beta2
    float pk_unpaired_penalty;     // beta3
    float pk_band_penalty;
    float pk_stack_span;
    float pk_interior_span;
    float multiloop_penalty_pk;
    float multiloop_paired_penalty_pk;
    float multiloop_unpaired_penalty_pk;
    float max_asymmetry;
    float salt_correction;
    float loop_greater30;
    float hairpin_GGG;
    float intermolecular_initiation;
} EnergyParms;

enum base_t { BASE_N = 0, BASE_A, BASE_C, BASE_G, BASE_U };
enum pair_t { PAIR_AU = 0, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG, PAIR_OTHER = -1 };

typedef Eigen::Tensor<float, 4> tensorN4;

class RNA
{
    public:
    RNA(string name, string seq);

    float  get_pij(int i, int j);
    string get_seq(void) const;
    uint   get_RNA_length(void) const;
    void   print_basepair_p_matrix(float theta) const;
    bool   allowed_basepair(size_t u, size_t v) const;
    bool   is_wc_basepair(size_t u, size_t v) const;

    private:
    base_t                                   base_type(char x) const;
    pair_t                                   pair_type(int i, int j) const;
    pair_t                                   pair_type(int i) const;    // assuming Watson-Crick pair
    void                                     load_default_parameters(void);
    float                                    GIL_asymmetry(uint l1, uint l2) const;
    float                                    Gpenalty(uint i, uint j) const;
    float                                    GIL(uint i, uint d, uint e, uint j, bool pk) const;
    float                                    GHL(uint i, uint j) const;
    float                                    GIL_mismatch(uint i, uint j, uint k, uint l) const;
    float                                    GIL_mismatch(uint i, uint j) const;
    float                                    Gloop(uint l) const;
    vector<MatrixXf>                         compute_partition_function_noPK_ON4(void);    // McCaskill
    vector<MatrixXf>                         compute_partition_function_noPK_ON3(void);    // McCaskill + fastILloops
    pair<vector<MatrixXf>, vector<tensorN4>> compute_partition_function_PK_ON8(void);      // McCaskill + PK
    pair<vector<MatrixXf>, vector<tensorN4>> compute_partition_function_PK_ON5(void);    // McCaskill + PK + fastILloops
    MatrixXf                                 compute_posterior_noPK_ON4(bool fast);      // probas
    MatrixXf                                 compute_posterior_PK_ON6(bool fast);        // probas + PK
    void compute_basepair_probabilities(vector<float> bp_proba, vector<int> offset, bool pk, bool fast);

    multi_array<pair_t, 2> pair_map;
    EnergyParms            nrjp_;    // energy parameters loaded from file or rna1995.h
    string                 name_;    // name of the rna
    string                 seq_;     // sequence of the rna with chars
    vector<base_t>         bseq_;    // sequence of the rna with base_ts
    uint                   n_;       // length of the rna
    vector<vector<float>>  pij_;     // vector of probabilities
};

inline float  RNA::get_pij(int i, int j) { return pij_[i][j]; }
inline uint   RNA::get_RNA_length() const { return n_; }
inline pair_t RNA::pair_type(int i, int j) const { return pair_map[bseq_[i]][bseq_[j]]; }
inline string RNA::get_seq(void) const { return seq_; }
inline bool   RNA::is_wc_basepair(size_t u, size_t v) const { return pair_type(u, v) != PAIR_OTHER; }

#endif
