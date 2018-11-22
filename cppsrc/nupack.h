/*
 * $Id$
 *
 * Copyright (C) 2010 Kengo Sato
 *
 * This file comes from IPknot.
 *
 */

#include <boost/multi_array.hpp>
#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::vector;

#define kB 0.00198717              // Boltzmann constant in kcal/mol/K
#define ZERO_C_IN_KELVIN 273.15    // Zero degrees C in Kelvin
#define AVOGADRO 6.022e23          // Avogadro's number

typedef float energy_t;

class DPtable2
{
  public:
    DPtable2() : V_(), N_(0) {}
    void resize(int n)
    {
        N_ = n;
        V_.resize(N_ * (N_ + 1) / 2 + (N_ + 1));
    }
    void          fill(const float& v) { std::fill(V_.begin(), V_.end(), v); }
    float&       operator()(int i, int j) { return V_[index(i, j)]; }
    const float& operator()(int i, int j) const { return V_[index(i, j)]; }

  private:
    int index(int i, int j) const
    {
        assert(j <= N_);
        return j == i - 1 ? N_ * (N_ + 1) / 2 + i : i * N_ + j - i * (1 + i) / 2;
    }
    std::vector<float> V_;
    int                 N_;
};

class DPtable4
{
  public:
    DPtable4() : V_(), N_(0) {}
    void resize(int n)
    {
        N_ = n;
        std::cout << V_.max_size() << " - " << N_ << " - " <<  sizeof(float) * static_cast<unsigned long>(N_) * (N_ - 1) * (N_ - 2) * (N_ - 3) / 2 / 3 / 4 << std::endl;
        V_.resize(static_cast<unsigned long>(N_) * (N_ - 1) * (N_ - 2) * (N_ - 3) / 2 / 3 / 4); // This number can be HUGE
        std::cout << "c'est toi qui bad_allocque ?" << std::endl;         
    }
    void          fill(const float& v) { std::fill(V_.begin(), V_.end(), v); }
    float&       operator()(int i, int d, int e, int j) { return V_[index(i, d, e, j)]; }
    const float& operator()(int i, int d, int e, int j) const { return V_[index(i, d, e, j)]; }

  private:
    int index(int h, int r, int m, int s) const
    {
        int n  = N_;
        int h2 = h * h;
        int h3 = h2 * h;
        int h4 = h3 * h;
        int m2 = m * m;
        int n2 = n * n;
        int n3 = n2 * n;
        int r2 = r * r;
        int r3 = r2 * r;
        assert(h <= r);
        assert(r <= m);
        assert(m <= s);
        assert(s <= N_);
        return (h == r && m == s) ? V_.size() - 1 :
                                    (
                                    -24 - 50 * h - 35 * h2 - 10 * h3 - h4 - 36 * m - 12 * m2 + 12 * n + 70 * h * n +
                                    30 * h2 * n + 4 * h3 * n + 24 * m * n - 12 * n2 - 30 * h * n2 - 6 * h2 * n2 + 4 * h * n3 +
                                    44 * r - 48 * n * r + 12 * n2 * r + 24 * r2 - 12 * n * r2 + 4 * r3 + 24 * s) /
                                    24;
    }
    std::vector<float> V_;
    int                 N_;
};

class DPtableX
{
  public:
    DPtableX() : V_(), N_(0), D_(0) {}
    void resize(int d, int n)
    {
        N_         = n;
        D_         = d;
        int max_sz = 0;
        for (int i = d; i < d + 3; ++i) max_sz = std::max(max_sz, (N_ - i) * (i - 5) * (i - 1) * (i - 2) / 2);
        V_.resize(max_sz);
    }
    void          fill(const float& v) { std::fill(V_.begin(), V_.end(), v); }
    float&       operator()(int i, int d, int e, int s) { return V_[index(i, d, e, s)]; }
    const float& operator()(int i, int d, int e, int s) const { return V_[index(i, d, e, s)]; }
    void          swap(DPtableX& x)
    {
        std::swap(V_, x.V_);
        std::swap(N_, x.N_);
        std::swap(D_, x.D_);
    }

  private:
    int index(int i, int h1, int m1, int s) const
    {
        int d      = D_;
        int d1d2   = (d - 1) * (d - 2);
        int d5     = d - 5;
        int h1_i_1 = h1 - i - 1;
        assert(i + d < N_);
        assert(d - 6 >= s);
        assert(i < h1);
        return i * d5 * d1d2 / 2 + s * d1d2 / 2 + h1_i_1 * (d - 1) - h1_i_1 * (h1 - i) / 2 + m1 - h1 - 1;
    }

    std::vector<float> V_;
    int                 N_;
    int                 D_;
};

class Nupack
{
  public:
    Nupack();
    void        load_sequence(const string& s);
    void        load_parameters_fm363(const vector<float>& v);
    void        load_default_parameters(/*int which*/);
    bool        load_parameters(const char* filename);
    void        dump_parameters(std::ostream& os) const;
    float calculate_partition_function();
    void        calculate_posterior();
    void        get_posterior(vector<float>& bp, vector<int>& offset) const;
    void        get_posterior(vector<float>& bp1, vector<float>& bp2, vector<int>& offset) const;

  private:
    void fastiloops(int i, int j, DPtable4& Qg, DPtableX& Qx, DPtableX& Qx2);
    void fastiloops_pr(int i, int j, DPtable4& Qg, DPtableX& Qx, DPtableX& Qx2, DPtable4& Pg, DPtableX& Px, DPtableX& Px2);
    energy_t score_hairpin(int i, int j) const;
    energy_t score_loop(int l) const;
    energy_t score_interior(int i, int d, int e, int j, bool pk) const;
    energy_t score_interior_mismatch(int i, int j) const;
    energy_t score_interior_mismatch(int i, int j, int k, int l) const;
    energy_t score_interior_asymmetry(int l1, int l2) const;
    energy_t score_multiloop(bool pk) const;
    energy_t score_multiloop_paired(int n, bool pk) const;
    energy_t score_multiloop_unpaired(int n, bool pk) const;
    energy_t score_at_penalty(int i, int j) const;
    energy_t score_dangle(int i, int j) const;
    energy_t score_pk() const;
    energy_t score_pk_multiloop() const;
    energy_t score_pk_pk() const;
    energy_t score_pk_paired(int n) const;
    energy_t score_pk_unpaired(int n) const;
    energy_t score_pk_band(int n) const;
    int      base(char x) const;
    bool     allow_paired(int i, int j) const;
    bool     wc_pair(int i, int j) const;
    int      pair_type(int i, int j) const;
    int      pair_type(int i) const;

    vector<int>                base_map;
    boost::multi_array<int, 2> pair_map;
    vector<int>                seq;
    int                        N;
    float                      RT;
    DPtable2                   Q;
    DPtable2                   Qb;
    DPtable2                   Qm;
    DPtable2                   Qp;
    DPtable2                   Qz;
    DPtable4                   Qg;
    DPtable4                   Qgl;
    DPtable4                   Qgr;
    DPtable4                   Qgls;
    DPtable4                   Qgrs;
    DPtable2                   P;
    DPtable2                   Pb;
    DPtable2                   Pm;
    DPtable2                   Pp;
    DPtable2                   Pz;
    DPtable2                   Pbg;
    DPtable4                   Pg;
    DPtable4                   Pgl;
    DPtable4                   Pgr;
    DPtable4                   Pgls;
    DPtable4                   Pgrs;

    // energy parameters
    energy_t hairpin37[30];
    energy_t bulge37[30];
    energy_t interior37[30];
    energy_t stack37[6][6];
    energy_t int11_37[6][6][4][4];
    energy_t int21_37[6][4][4][6][4];
    energy_t int22_37[6][6][4][4][4][4];
    energy_t dangle3_37[6][4];
    energy_t dangle5_37[6][4];
    energy_t triloop37[4][4][4][4][4];
    energy_t tloop37[4][4][4][4][4][4];
    energy_t mismatch_hairpin37[4][4][6];
    energy_t mismatch_interior37[4][4][6];
    energy_t asymmetry_penalty[4];
    energy_t polyC_penalty, polyC_slope, polyC_int;
    energy_t at_penalty;
    energy_t multiloop_penalty;             // alpha1
    energy_t multiloop_paired_penalty;      // alpha2
    energy_t multiloop_unpaired_penalty;    // alpha3
    energy_t pk_penalty;                    // beta1
    energy_t pk_multiloop_penalty;          // beta1m
    energy_t pk_pk_penalty;                 // beta1p
    energy_t pk_paired_penalty;             // beta2
    energy_t pk_unpaired_penalty;           // beta3
    energy_t pk_band_penalty;
    energy_t pk_stack_span;
    energy_t pk_interior_span;
    energy_t multiloop_penalty_pk;
    energy_t multiloop_paired_penalty_pk;
    energy_t multiloop_unpaired_penalty_pk;
    energy_t max_asymmetry;
    energy_t salt_correction;
    energy_t loop_greater30;
    energy_t hairpin_GGG;
    float    intermolecular_initiation;
};
