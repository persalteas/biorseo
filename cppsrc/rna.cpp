#include "rna.h"
#include "rna1995.h"
#include <Eigen/Core>
#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace Eigen;

using std::abs;
using std::cerr;
using std::cout;
using std::endl;

RNA::RNA(string name, string seq)
{
    name_ = name;
    seq_  = seq;
    n_    = seq_.size();
    bseq_.reserve(n_);
    vector<char> unknown_chars;
    bool         contains_T = false;
    for (char c : seq) {
        if (c == 'T' or c == 't') {
            c          = 'U';
            contains_T = true;
        }
        if (base_type(c) == BASE_N) {
            unknown_chars.push_back(c);
        }
        bseq_.push_back(base_type(c));
    }
    if (contains_T) cout << "\tWARNING: Thymines automatically replaced by uraciles.";
    if (unknown_chars.size() > 0) {
        cout << "\tWARNING: Unknown chars in input sequence ignored : ";
        for (char c : unknown_chars) cout << c << " ";
        cout << endl;
    }
    cout << "\t>Sequence formatted" << endl;

    /*define pij_*/
    vector<float> bp_proba;
    vector<int>   offset;
    // load_parameters("rna1999.dG"); // to load custom parameters
    load_default_parameters();
    cout << "\t>default parameters loaded (Serra and Turner, 1995)" << endl;
    cout << "\t>computing pairing probabilities..." << endl;
    try {    // Depending on the input RNA size, the RAM amount might too large to handle for your poor laptop
        compute_partition_function();
    } catch (std::exception& e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    compute_posterior();
    get_posterior(bp_proba, offset);

    pij_ = vector<vector<float>>(n_, vector<float>(n_));
    for (int i = 1; i <= n_; i++) {
        for (int j = 1; j <= n_; j++) {
            pij_[i - 1][j - 1] = bp_proba[offset[i] + j];
        }
    }
    cout << "\t>pairing probabilities defined" << endl;
}

void RNA::print_basepair_p_matrix(float theta) const
{
    cout << endl << endl;
    cout << "\t=== -log10(p(i,j)) for each pair (i,j) of nucleotides: ===" << endl << endl;
    cout << "\t " << seq_ << endl;
    uint i = 0;
    for (auto u : pij_) {
        cout << "\t" << seq_[i];
        for (auto v : u) {
            if (v < 5e-10)
                cout << " ";
            else if (v > theta)
                cout << "\033[0;32m" << int(-log10(v)) << "\033[0m";
            else
                cout << int(-log10(v));
        }
        cout << endl;
        i++;
    }
    cout << endl << "\t\033[0;32mgreen\033[0m basepairs are kept as decision variables." << endl << endl;
}

void RNA::load_default_parameters()
{
    int p[] = {PAIR_AU, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG};
    int b[] = {BASE_A - 1, BASE_C - 1, BASE_G - 1, BASE_U - 1};

    const int* v = &thermo_params[0];

    // stack
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j) nrjp_.stack37[p[i]][p[j]] = *(v++) / 100.0;

    // hairpin
    for (int i = 0; i < 30; ++i) nrjp_.hairpin37[i] = *(v++) / 100.0;

    // bulge
    for (int i = 0; i < 30; ++i) nrjp_.bulge37[i] = *(v++) / 100.0;

    // interior
    for (int i = 0; i < 30; ++i) nrjp_.interior37[i] = *(v++) / 100.0;

    // asymmetry panelties
    for (int i = 0; i < 4; ++i) nrjp_.asymmetry_penalty[i] = *(v++) / 100.0;

    // mismatch hairpin
    for (int i = 0; i != 4; ++i)
        for (int j = 0; j != 4; ++j)
            for (int k = 0; k != 6; ++k) nrjp_.mismatch_hairpin37[b[i]][b[j]][p[k]] = *(v++) / 100.0;

    // mismatch interior
    for (int i = 0; i != 4; ++i)
        for (int j = 0; j != 4; ++j)
            for (int k = 0; k != 6; ++k) nrjp_.mismatch_interior37[b[i]][b[j]][p[k]] = *(v++) / 100.0;

    // dangle5
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 4; ++j) nrjp_.dangle5_37[p[i]][b[j]] = *(v++) / 100.0;

    // dangle3
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 4; ++j) nrjp_.dangle3_37[p[i]][b[j]] = *(v++) / 100.0;

    // multiloop penalties
    nrjp_.a1 = *(v++) / 100.0;
    nrjp_.a2 = *(v++) / 100.0;
    nrjp_.a3 = *(v++) / 100.0;

    // AT terminate penalties
    nrjp_.at_penalty = *(v++) / 100.0;

    // interior loops 1x1
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int k = 0; k != 4; ++k)
                for (int l = 0; l != 4; ++l) nrjp_.int11_37[p[i]][p[j]][b[k]][b[l]] = *(v++) / 100.0;

    // interior loops 2x2
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int m = 0; m != 4; ++m)
                for (int n = 0; n != 4; ++n)
                    for (int k = 0; k != 4; ++k)
                        for (int l                                             = 0; l != 4; ++l)
                            nrjp_.int22_37[p[i]][p[j]][b[m]][b[l]][b[n]][b[k]] = *(v++) / 100.0;

    // interior loops 1x2
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int m = 0; m != 4; ++m)
                for (int k = 0; k != 4; ++k)
                    for (int l = 0; l != 4; ++l) nrjp_.int21_37[p[i]][b[k]][b[m]][p[j]][b[l]] = *(v++) / 100.0;

    // nrjp_.polyC hairpin parameters
    nrjp_.polyC_penalty = *(v++) / 100.0;
    nrjp_.polyC_slope   = *(v++) / 100.0;
    nrjp_.polyC_int     = *(v++) / 100.0;

    // pseudoknot energy parameters
    nrjp_.pk_penalty                    = *(v++) / 100.0;
    nrjp_.pk_paired_penalty             = *(v++) / 100.0;
    nrjp_.pk_unpaired_penalty           = *(v++) / 100.0;
    nrjp_.pk_multiloop_penalty          = *(v++) / 100.0;
    nrjp_.pk_pk_penalty                 = *(v++) / 100.0;
    nrjp_.pk_band_penalty               = 0.0;
    nrjp_.pk_stack_span                 = 1.0;
    nrjp_.pk_interior_span              = 1.0;
    nrjp_.multiloop_penalty_pk          = nrjp_.a1;
    nrjp_.multiloop_paired_penalty_pk   = nrjp_.a2;
    nrjp_.multiloop_unpaired_penalty_pk = nrjp_.a3;

    // BIMOLECULAR TERM
    nrjp_.intermolecular_initiation = *(v++) / 100.0;

    // triloops
    // std::fill(nrjp_.triloop37.data(), nrjp_.triloop37.data()+nrjp_.triloop37.num_elements(),
    // 0.0);
    std::fill(&nrjp_.triloop37[0][0][0][0][0], &nrjp_.triloop37[0][0][0][0][0] + 4 * 4 * 4 * 4 * 4, 0.0);
    for (int i = 0; triloops[i].s != NULL; ++i) {
        int         v    = triloops[i].e;
        const char* loop = triloops[i].s;
        vector<int> idx(5);
        for (int i = 0; i != 5; ++i) idx[i] = base_type(loop[i]) - 1;
        // nrjp_.triloop37(idx) = v/100.0;
        nrjp_.triloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]] = v / 100.0;
    }

    // tloops
    // std::fill(nrjp_.tloop37.data(), nrjp_.tloop37.data()+nrjp_.tloop37.num_elements(), 0.0);
    std::fill(&nrjp_.tloop37[0][0][0][0][0][0], &nrjp_.tloop37[0][0][0][0][0][0] + 4 * 4 * 4 * 4 * 4 * 4, 0.0);
    for (int i = 0; tetraloops[i].s != NULL; ++i) {
        int         v    = tetraloops[i].e;
        const char* loop = tetraloops[i].s;
        vector<int> idx(6);
        for (int i = 0; i != 6; ++i) idx[i] = base_type(loop[i]) - 1;
        // nrjp_.tloop37(idx) = v/100.0;
        nrjp_.tloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]] = v / 100.0;
    }
}

float RNA::score_interior_mismatch(int i, int j, int k, int l) const
{
    return nrjp_.mismatch_interior37[seq_[k] - 1][seq_[l] - 1][pair_type(i, j)];
}

float RNA::score_interior_mismatch(int i, int j) const
{
    return nrjp_.mismatch_interior37[BASE_N][BASE_N][pair_type(i, j)];
}


float RNA::score_at_penalty(int i, int j) const
{
    return pair_type(i, j) == PAIR_AU || pair_type(i, j) == PAIR_UA ? nrjp_.at_penalty : 0;
}


float_t RNA::score_interior_asymmetry(int l1, int l2) const
{
    float e         = 0.0;
    int   size      = l1 + l2;
    int   asymmetry = abs(l1 - l2);
    e += size <= 30 ? nrjp_.interior37[size - 1] : nrjp_.interior37[30 - 1] + nrjp_.loop_greater30 * log(size / 30.0);

    // asymmetry penalty
    e += std::min(nrjp_.max_asymmetry, asymmetry * nrjp_.asymmetry_penalty[std::min(4, std::min(l1, l2)) - 1]);

    return e;
}

base_t RNA::base_type(char x) const
{
    if (x == 'a' or x == 'A') return BASE_A;
    if (x == 'c' or x == 'C') return BASE_A;
    if (x == 'g' or x == 'G') return BASE_A;
    if (x == 'u' or x == 'U') return BASE_A;
    return BASE_N;
}

pair_t RNA::pair_type(int i) const
{
    // assume Watson-Crick pairs
    switch (bseq_[i]) {
    case BASE_A: return PAIR_AU; break;
    case BASE_C: return PAIR_CG; break;
    case BASE_G: return PAIR_GC; break;
    case BASE_U: return PAIR_UA;
    }
}

float RNA::Ghairpin(uint i, uint j) const
{
    float e     = 0.0;
    bool  polyC = true;
    for (int k = i + 1; k < j; ++k) {
        if (seq_[k] != BASE_C) {
            polyC = false;
            break;
        }
    }

    int size = j - i - 1;

    assert(size >= 3);
    // assert(allow_basepair(i, j));

    e += size <= 30 ? nrjp_.hairpin37[size - 1] : nrjp_.hairpin37[30 - 1] + nrjp_.loop_greater30 * log(size / 30.0);

    if (size == 3) {
        e += score_at_penalty(i, j);
        e += nrjp_.triloop37[seq_[i] - 1][seq_[i + 1] - 1][seq_[i + 2] - 1][seq_[j - 1] - 1][seq_[j] - 1];
        if (polyC) e += nrjp_.polyC_penalty;
        if (seq_[i + 1] == BASE_G && seq_[i + 2] == BASE_G && seq_[j - 1] == BASE_G) e += nrjp_.hairpin_GGG;
    } else if (size == 4) {
        e += nrjp_.tloop37[seq_[i] - 1][seq_[i + 1] - 1][seq_[i + 2] - 1][seq_[j - 2] - 1][seq_[j - 1] - 1][seq_[j] - 1];
        e += nrjp_.mismatch_hairpin37[seq_[i + 1] - 1][seq_[j - 1] - 1][pair_type(i, j)];
        if (polyC) e += nrjp_.polyC_slope * size + nrjp_.polyC_int;
    } else /*if (size>4)*/
    {
        e += nrjp_.mismatch_hairpin37[seq_[i + 1] - 1][seq_[j - 1] - 1][pair_type(i, j)];
        if (polyC) e += nrjp_.polyC_slope * size + nrjp_.polyC_int;
    }
    return e;
}


float RNA::Ginterior(uint i, uint h, uint m, uint j, bool pk) const
{
    int   l1   = h - i - 1;
    int   l2   = j - m - 1;
    int   size = l1 + l2;
    float e    = 0;

    // helix
    if (size == 0) {
        return nrjp_.stack37[pair_type(i, j)][pair_type(h, m)] * (pk ? nrjp_.pk_stack_span : 1.0);
    }

    // bulge
    else if (l1 == 0 || l2 == 0) {
        e += size <= 30 ? nrjp_.bulge37[size - 1] : nrjp_.bulge37[30 - 1] + nrjp_.loop_greater30 * log(size / 30.0);

        if (l1 + l2 == 1)    // single bulge...treat as a stacked region
        {
            e += nrjp_.stack37[pair_type(i, j)][pair_type(h, m)];
            e -= nrjp_.salt_correction;
        } else {
            e += score_at_penalty(i, j);
            e += score_at_penalty(h, m);
        }
    }

    // interior loop
    else if (l1 > 0 && l2 > 0) {
        int asymmetry = abs(l1 - l2);
        if (asymmetry > 1 || size > 4) {
            e += score_interior_asymmetry(l1, l2);
            if (l1 > 1 && l2 > 1) {
                e += score_interior_mismatch(m, h, m + 1, h - 1);
                e += score_interior_mismatch(i, j, i + 1, j - 1);
            } else if (l1 == 1 || l2 == 1) {
                e += score_interior_mismatch(m, h);
                e += score_interior_mismatch(i, j);
            } else {
                assert(!"unclassified interior loop");
                exit(1);
            }
        } else if (l1 == 1 && l2 == 1)
            e += nrjp_.int11_37[pair_type(i, j)][pair_type(h, m)][seq_[i + 1] - 1][seq_[j - 1] - 1];
        else if (l1 == 2 && l2 == 2)
            e += nrjp_.int22_37[pair_type(i, j)][pair_type(h, m)][seq_[i + 1] - 1][seq_[j - 1] - 1][seq_[i + 2] - 1]
                               [seq_[j - 2] - 1];
        else if (l1 == 1 && l2 == 2)
            e += nrjp_.int21_37[pair_type(i, j)][seq_[j - 2] - 1][seq_[i + 1] - 1][pair_type(h, m)][seq_[j - 1] - 1];
        else if (l1 == 2 && l2 == 1)
            e += nrjp_.int21_37[pair_type(m, h)][seq_[i + 1] - 1][seq_[j - 1] - 1][pair_type(j, i)][seq_[i + 2] - 1];
        else {
            assert(!"error in tabulated interior loop");
            exit(EXIT_FAILURE);
        }
    } else {
        assert(!"improperly classifed interior loop");
        exit(EXIT_FAILURE);
    }
    return e * (pk ? nrjp_.pk_interior_span : 1.0);
}

float RNA::compute_partition_function(void)
{
    // This is the o(N⁴) algorithm from Dirks & Pierce, 2003
    // Basically is the same than McCaskill, 1990
    // Gmultiloop is approximated by a1 + k*a2 + u*a3 (a,b,c in McCaskill)
    // Pseudoknots supposed impossible

    float RT = kB * AVOGADRO * (ZERO_C_IN_KELVIN + 37.0);
    float a1 = nrjp_.a1;
    float a2 = nrjp_.a2;
    float a3 = nrjp_.a3;

    // O(3N²) space
    MatrixXf Q  = MatrixXf::Zero(n_, n_);
    MatrixXf Qb = MatrixXf::Zero(n_, n_);
    MatrixXf Qm = MatrixXf::Zero(n_, n_);
    for (uint i = 1; i < n_; i++) Q(i, i - 1) = 1.0;

    for (uint l = 1; l <= n_; l++)    // Consider subsequences of growing sizes until n_
    {
        // This loop on i should be parallelized
        for (uint i = 0; i < n_ - l + 1; i++) {
            uint j = i + l - 1;    // Consider the subsequence [i,j] of length l

            // Qb recursion
            Qb(i, j) = exp(-Ghairpin(i, j) / RT);
            for (uint d = i + 1; d <= j - 5; d++)    // loop over all possible rightmost basepairs (d,e)
            {
                for (uint e = d + 4; e <= j - 1; e++) {
                    Qb(i, j) += exp(-Ginterior(i, d, e, j, false) / RT) * Qb(d, e);
                    Qb(i, j) += Qm(i + 1, d - 1) * Qb(d, e) * exp(-(a1 + 2 * a2 + (j - e - 1) * a3) / RT);
                }
            }

            // Q and Qm recursion
            Q(i, j) = 1.0;                         // if empty (no basepairs between i and j)
            for (uint d = i; d <= j - 4; d++) {    // loop over all possible rightmost basepairs (d,e)
                for (uint e = d + 4; e <= j; e++) {
                    Q(i, j) += Q(i, d - 1) * Qb(d, e);
                    Qm(i, j) += exp(-(a2 + a3 * (d - i) + a3 * (j - e)) / RT) * Qb(d, e);
                    Qm(i, j) += Qm(i, d - 1) * Qb(d, e) * exp(-(a2 + a3 * (j - e)) / RT);
                }
            }
        }
    }
    return Q(0, n_ - 1);    // Partition function is Q(1,N)
}