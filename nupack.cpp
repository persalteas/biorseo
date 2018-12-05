#include "nupack.h"
#include "rna1995.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

enum { A = BASE_A - 1, C = BASE_C - 1, G = BASE_G - 1, U = BASE_U - 1 };
enum { AU = PAIR_AU, CG = PAIR_CG, GC = PAIR_GC, UA = PAIR_UA, GU = PAIR_GU, UG = PAIR_UG };

#define EXP expl
#define LOG logl

Nupack::Nupack()
: base_map('z' - 'a' + 1), pair_map(boost::extents[5][5]), RT(kB * (ZERO_C_IN_KELVIN + 37)), salt_correction(0),
  loop_greater30(1.079 /*=1.75*RT*/), hairpin_GGG(0.0)
{
    std::fill(base_map.begin(), base_map.end(), (int)BASE_N);
    base_map['a' - 'a'] = BASE_A;
    base_map['c' - 'a'] = BASE_C;
    base_map['g' - 'a'] = BASE_G;
    base_map['u' - 'a'] = BASE_U;
    base_map['t' - 'a'] = BASE_U;

    std::fill(pair_map.data(), pair_map.data() + pair_map.num_elements(), -1);
    pair_map[BASE_A][BASE_U] = PAIR_AU;
    pair_map[BASE_U][BASE_A] = PAIR_UA;
    pair_map[BASE_C][BASE_G] = PAIR_CG;
    pair_map[BASE_G][BASE_C] = PAIR_GC;
    pair_map[BASE_G][BASE_U] = PAIR_GU;
    pair_map[BASE_U][BASE_G] = PAIR_UG;
}

int Nupack::base(char x) const
{
    if (x >= 'a' && x <= 'z')
        return base_map[x - 'a'];
    else if (x >= 'A' && x <= 'Z')
        return base_map[x - 'A'];
    else
        return BASE_N;
}

int Nupack::pair_type(int i, int j) const { return pair_map[seq[i]][seq[j]]; }

int Nupack::pair_type(int i) const
{
    // assume Watson-Crick pairs
    switch (seq[i]) {
    case BASE_A: return PAIR_AU; break;
    case BASE_C: return PAIR_CG; break;
    case BASE_G: return PAIR_GC; break;
    case BASE_U: return PAIR_UA; break;
    }
    return -1;
}

bool Nupack::wc_pair(int i, int j) const { return pair_type(i, j) != PAIR_GU && pair_type(i, j) != PAIR_UG; }

bool Nupack::allow_paired(int i, int j) const { return j - i - 1 >= 3 && pair_type(i, j) >= 0; }

void Nupack::load_sequence(const string& s)
{
    N = s.size();
    seq.resize(N);
    for (int i = 0; i != N; ++i) seq[i] = base(s[i]);
}

bool Nupack::load_parameters(const char* file)
{
    int p[] = {PAIR_AU, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG};
    int b[] = {BASE_A - 1, BASE_C - 1, BASE_G - 1, BASE_U - 1};

    std::ifstream is(file);
    if (!is) return false;

    string line;

    // stack
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 6; ++i) {
        std::istringstream ss(line);
        for (int j = 0; j != 6; ++j) {
            int v;
            ss >> v;
            stack37[p[i]][p[j]] = v / 100.0;
        }
        std::getline(is, line);
    }

    // hairpin
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    {
        int                i, j, v;
        std::istringstream ss(line);
        for (i = 0; i < 30 && ss; ++i) {
            ss >> v;
            hairpin37[i] = v / 100.0;
        }
        for (j = i; j < 30; ++j) hairpin37[j] = hairpin37[i - 1] + loop_greater30 * LOG((j + 1) / (1.0 * i));
    }

    // bulge
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    {
        int                i, j, v;
        std::istringstream ss(line);
        for (i = 0; i < 30 && ss; ++i) {
            ss >> v;
            bulge37[i] = v / 100.0;
        }
        for (j = i; j < 30; ++j) bulge37[j] = bulge37[i - 1] + loop_greater30 * LOG((j + 1) / (1.0 * i));
    }

    // interior
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    {
        int                i, j, v;
        std::istringstream ss(line);
        for (i = 0; i < 30 && ss; ++i) {
            ss >> v;
            interior37[i] = v / 100.0;
        }
        for (j = i; j < 30; ++j) interior37[j] = interior37[i - 1] + loop_greater30 * LOG((j + 1) / (1.0 * i));
    }

    // asymmetry panelties
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    {
        int                v;
        std::istringstream ss(line);
        for (int i = 0; i < 4; ++i) {
            ss >> v;
            asymmetry_penalty[i] = v / 100.0;
        }
        ss >> v;
        max_asymmetry = v / 100.0;
    }

    // triloops
    // std::fill(triloop37.data(), triloop37.data()+triloop37.num_elements(),
    // 0.0);
    std::fill(&triloop37[0][0][0][0][0], &triloop37[0][0][0][0][0] + 4 * 4 * 4 * 4 * 4, 0.0);
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    while (line[0] != '>') {
        int  v;
        char loop[256];
        sscanf(line.c_str(), "%s %d", loop, &v);
        vector<int> idx(5);
        for (int i = 0; i != 5; ++i) idx[i] = base(loop[i]) - 1;
        // triloop37(idx) = v/100.0;
        triloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]] = v / 100.0;
        std::getline(is, line);
    }

    // tloops
    // std::fill(tloop37.data(), tloop37.data()+tloop37.num_elements(), 0.0);
    std::fill(&tloop37[0][0][0][0][0][0], &tloop37[0][0][0][0][0][0] + 4 * 4 * 4 * 4 * 4 * 4, 0.0);
    std::getline(is, line);
    while (line[0] == '>') std::getline(is, line);
    while (line[0] != '>') {
        int  v;
        char loop[256];
        sscanf(line.c_str(), "%s %d", loop, &v);
        vector<int> idx(6);
        for (int i = 0; i != 6; ++i) idx[i] = base(loop[i]) - 1;
        // tloop37(idx) = v/100.0;
        tloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]] = v / 100.0;
        std::getline(is, line);
    }

    // mismatch hairpin
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 4; ++i) {
        for (int j = 0; j != 4; ++j) {
            std::istringstream ss(line);
            for (int k = 0; k != 6; ++k) {
                int v;
                ss >> v;
                mismatch_hairpin37[b[i]][b[j]][p[k]] = v / 100.0;
            }
            std::getline(is, line);
        }
    }

    // mismatch interior
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 4; ++i) {
        for (int j = 0; j != 4; ++j) {
            std::istringstream ss(line);
            for (int k = 0; k != 6; ++k) {
                int v;
                ss >> v;
                mismatch_interior37[b[i]][b[j]][p[k]] = v / 100.0;
            }
            std::getline(is, line);
        }
    }

    // dangle5
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 6; ++i) {
        std::istringstream ss(line);
        for (int j = 0; j != 4; ++j) {
            int v;
            ss >> v;
            dangle5_37[p[i]][b[j]] = v / 100.0;
        }
        std::getline(is, line);
    }

    // dangle3
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 6; ++i) {
        std::istringstream ss(line);
        for (int j = 0; j != 4; ++j) {
            int v;
            ss >> v;
            dangle3_37[p[i]][b[j]] = v / 100.0;
        }
        std::getline(is, line);
    }

    // multiloop penalties
    while (line[0] == '>') std::getline(is, line);
    {
        int                v;
        std::istringstream ss(line);
        ss >> v;
        multiloop_penalty = v / 100.0;
        ss >> v;
        multiloop_paired_penalty = v / 100.0;
        ss >> v;
        multiloop_unpaired_penalty = v / 100.0;
        std::getline(is, line);
    }

    // AT terminate penalties
    while (line[0] == '>') std::getline(is, line);
    {
        int                v;
        std::istringstream ss(line);
        ss >> v;
        at_penalty = v / 100.0;
        std::getline(is, line);
    }

    // interior loops 1x1
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 6; ++i) {
        for (int j = 0; j != 6; ++j) {
            std::getline(is, line);    // header
            for (int k = 0; k != 4; ++k) {
                std::istringstream ss(line);
                for (int l = 0; l != 4; ++l) {
                    int v;
                    ss >> v;
                    int11_37[p[i]][p[j]][b[k]][b[l]] = v / 100.0;
                }
                std::getline(is, line);
            }
        }
    }

    // interior loops 2x2
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 6; ++i) {
        for (int j = 0; j != 6; ++j) {
            for (int m = 0; m != 4; ++m) {
                for (int n = 0; n != 4; ++n) {
                    std::getline(is, line);    // header
                    for (int k = 0; k != 4; ++k) {
                        std::istringstream ss(line);
                        for (int l = 0; l != 4; ++l) {
                            int v;
                            ss >> v;
                            int22_37[p[i]][p[j]][b[m]][b[l]][b[n]][b[k]] = v / 100.0;
                        }
                        std::getline(is, line);
                    }
                }
            }
        }
    }

    // interior loops 1x2
    while (line[0] == '>') std::getline(is, line);
    for (int i = 0; i != 6; ++i) {
        for (int j = 0; j != 6; ++j) {
            for (int m = 0; m != 4; ++m) {
                std::getline(is, line);    // header
                for (int k = 0; k != 4; ++k) {
                    std::istringstream ss(line);
                    for (int l = 0; l != 4; ++l) {
                        int v;
                        ss >> v;
                        int21_37[p[i]][b[k]][b[m]][p[j]][b[l]] = v / 100.0;
                    }
                    std::getline(is, line);
                }
            }
        }
    }

    // polyC hairpin parameters
    while (line[0] == '>') std::getline(is, line);
    {
        int                v;
        std::istringstream ss(line);
        ss >> v;
        polyC_penalty = v / 100.0;
        ss >> v;
        polyC_slope = v / 100.0;
        ss >> v;
        polyC_int = v / 100.0;
        std::getline(is, line);
    }

    // pseudoknot energy parameters
    while (line[0] == '>') std::getline(is, line);
    {
        int                v;
        std::istringstream ss(line);
        ss >> v;
        pk_penalty = v / 100.0;
        ss >> v;
        pk_paired_penalty = v / 100.0;
        ss >> v;
        pk_unpaired_penalty = v / 100.0;
        ss >> v;
        pk_multiloop_penalty = v / 100.0;
        ss >> v;
        pk_pk_penalty = v / 100.0;
        std::getline(is, line);
    }
    pk_band_penalty               = 0.0;
    pk_stack_span                 = 1.0;
    pk_interior_span              = 1.0;
    multiloop_penalty_pk          = multiloop_penalty;
    multiloop_paired_penalty_pk   = multiloop_paired_penalty;
    multiloop_unpaired_penalty_pk = multiloop_unpaired_penalty;

    // BIMOLECULAR TERM
    while (line[0] == '>') std::getline(is, line);
    {
        int                v;
        std::istringstream ss(line);
        ss >> v;
        intermolecular_initiation = v / 100.0;
        std::getline(is, line);
    }

    return true;
}

void Nupack::load_default_parameters()
{
    int p[] = {PAIR_AU, PAIR_CG, PAIR_GC, PAIR_UA, PAIR_GU, PAIR_UG};
    int b[] = {BASE_A - 1, BASE_C - 1, BASE_G - 1, BASE_U - 1};

    const int* v = &thermo_params[0];

    // stack
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j) stack37[p[i]][p[j]] = *(v++) / 100.0;

    // hairpin
    for (int i = 0; i < 30; ++i) hairpin37[i] = *(v++) / 100.0;

    // bulge
    for (int i = 0; i < 30; ++i) bulge37[i] = *(v++) / 100.0;

    // interior
    for (int i = 0; i < 30; ++i) interior37[i] = *(v++) / 100.0;

    // asymmetry panelties
    for (int i = 0; i < 4; ++i) asymmetry_penalty[i] = *(v++) / 100.0;

    // mismatch hairpin
    for (int i = 0; i != 4; ++i)
        for (int j = 0; j != 4; ++j)
            for (int k = 0; k != 6; ++k) mismatch_hairpin37[b[i]][b[j]][p[k]] = *(v++) / 100.0;

    // mismatch interior
    for (int i = 0; i != 4; ++i)
        for (int j = 0; j != 4; ++j)
            for (int k = 0; k != 6; ++k) mismatch_interior37[b[i]][b[j]][p[k]] = *(v++) / 100.0;

    // dangle5
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 4; ++j) dangle5_37[p[i]][b[j]] = *(v++) / 100.0;

    // dangle3
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 4; ++j) dangle3_37[p[i]][b[j]] = *(v++) / 100.0;

    // multiloop penalties
    multiloop_penalty          = *(v++) / 100.0;
    multiloop_paired_penalty   = *(v++) / 100.0;
    multiloop_unpaired_penalty = *(v++) / 100.0;

    // AT terminate penalties
    at_penalty = *(v++) / 100.0;

    // interior loops 1x1
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int k = 0; k != 4; ++k)
                for (int l = 0; l != 4; ++l) int11_37[p[i]][p[j]][b[k]][b[l]] = *(v++) / 100.0;

    // interior loops 2x2
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int m = 0; m != 4; ++m)
                for (int n = 0; n != 4; ++n)
                    for (int k = 0; k != 4; ++k)
                        for (int l = 0; l != 4; ++l) int22_37[p[i]][p[j]][b[m]][b[l]][b[n]][b[k]] = *(v++) / 100.0;

    // interior loops 1x2
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int m = 0; m != 4; ++m)
                for (int k = 0; k != 4; ++k)
                    for (int l = 0; l != 4; ++l) int21_37[p[i]][b[k]][b[m]][p[j]][b[l]] = *(v++) / 100.0;

    // polyC hairpin parameters
    polyC_penalty = *(v++) / 100.0;
    polyC_slope   = *(v++) / 100.0;
    polyC_int     = *(v++) / 100.0;

    // pseudoknot energy parameters
    pk_penalty                    = *(v++) / 100.0;
    pk_paired_penalty             = *(v++) / 100.0;
    pk_unpaired_penalty           = *(v++) / 100.0;
    pk_multiloop_penalty          = *(v++) / 100.0;
    pk_pk_penalty                 = *(v++) / 100.0;
    pk_band_penalty               = 0.0;
    pk_stack_span                 = 1.0;
    pk_interior_span              = 1.0;
    multiloop_penalty_pk          = multiloop_penalty;
    multiloop_paired_penalty_pk   = multiloop_paired_penalty;
    multiloop_unpaired_penalty_pk = multiloop_unpaired_penalty;

    // BIMOLECULAR TERM
    intermolecular_initiation = *(v++) / 100.0;

    // triloops
    // std::fill(triloop37.data(), triloop37.data()+triloop37.num_elements(),
    // 0.0);
    std::fill(&triloop37[0][0][0][0][0], &triloop37[0][0][0][0][0] + 4 * 4 * 4 * 4 * 4, 0.0);
    for (int i = 0; triloops[i].s != NULL; ++i) {
        int         v    = triloops[i].e;
        const char* loop = triloops[i].s;
        vector<int> idx(5);
        for (int i = 0; i != 5; ++i) idx[i] = base(loop[i]) - 1;
        // triloop37(idx) = v/100.0;
        triloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]] = v / 100.0;
    }

    // tloops
    // std::fill(tloop37.data(), tloop37.data()+tloop37.num_elements(), 0.0);
    std::fill(&tloop37[0][0][0][0][0][0], &tloop37[0][0][0][0][0][0] + 4 * 4 * 4 * 4 * 4 * 4, 0.0);
    for (int i = 0; tetraloops[i].s != NULL; ++i) {
        int         v    = tetraloops[i].e;
        const char* loop = tetraloops[i].s;
        vector<int> idx(6);
        for (int i = 0; i != 6; ++i) idx[i] = base(loop[i]) - 1;
        // tloop37(idx) = v/100.0;
        tloop37[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]] = v / 100.0;
    }
}

void Nupack::dump_parameters(std::ostream& os) const
{
    // stack
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j) os << "stack37[" << i << "][" << j << "]=" << stack37[i][j] << std::endl;

    // hairpin
    for (int i = 0; i != 30; ++i) os << "hairpin37[" << i << "]=" << hairpin37[i] << std::endl;

    // bulge
    for (int i = 0; i != 30; ++i) os << "bulge37[" << i << "]=" << bulge37[i] << std::endl;

    // interior
    for (int i = 0; i != 30; ++i) os << "interior37[" << i << "]=" << interior37[i] << std::endl;

    // asymmetry
    for (int i = 0; i != 4; ++i) os << "asymmetry_penalty[" << i << "]=" << asymmetry_penalty[i] << std::endl;
    os << "max_asymmetry=" << max_asymmetry << std::endl;

    // triloops
    for (int i0 = 0; i0 != 4; ++i0)
        for (int i1 = 0; i1 != 4; ++i1)
            for (int i2 = 0; i2 != 4; ++i2)
                for (int i3 = 0; i3 != 4; ++i3)
                    for (int i4 = 0; i4 != 4; ++i4)
                        if (triloop37[i0][i1][i2][i3][i4] != 0.0)
                            os << "triloop37[" << i0 << "][" << i1 << "][" << i2 << "][" << i3 << "][" << i4
                               << "]=" << triloop37[i0][i1][i2][i3][i4] << std::endl;

    // tloops
    for (int i0 = 0; i0 != 4; ++i0)
        for (int i1 = 0; i1 != 4; ++i1)
            for (int i2 = 0; i2 != 4; ++i2)
                for (int i3 = 0; i3 != 4; ++i3)
                    for (int i4 = 0; i4 != 4; ++i4)
                        for (int i5 = 0; i5 != 4; ++i5)
                            if (tloop37[i0][i1][i2][i3][i4][i5] != 0.0)
                                os << "tloop37[" << i0 << "][" << i1 << "][" << i2 << "][" << i3 << "][" << i4 << "]["
                                   << i5 << "]=" << tloop37[i0][i1][i2][i3][i4][i5] << std::endl;

    // mismatch hairpin
    for (int i = 0; i != 4; ++i)
        for (int j = 0; j != 4; ++j)
            for (int k = 0; k != 6; ++k)
                os << "mismatch_hairpin37[" << i << "][" << j << "][" << k << "]=" << mismatch_hairpin37[i][j][k] << std::endl;

    // mismatch interior37
    for (int i = 0; i != 4; ++i)
        for (int j = 0; j != 4; ++j)
            for (int k = 0; k != 6; ++k)
                os << "mismatch_interior37[" << i << "][" << j << "][" << k << "]=" << mismatch_interior37[i][j][k] << std::endl;

    // dangle5
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 4; ++j) os << "dangle5_37[" << i << "][" << j << "]=" << dangle5_37[i][j] << std::endl;

    // dangle3
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 4; ++j) os << "dangle3_37[" << i << "][" << j << "]=" << dangle3_37[i][j] << std::endl;

    // multiloop penalties
    os << "multiloop_penalty=" << multiloop_penalty << std::endl
       << "multiloop_paired_penalty=" << multiloop_paired_penalty << std::endl
       << "multiloop_unpaired_penalty=" << multiloop_unpaired_penalty << std::endl;

    // AT terminate penalties
    os << "at_penalty=" << at_penalty << std::endl;

    // interior loops 1x1
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int k = 0; k != 4; ++k)
                for (int l = 0; l != 4; ++l)
                    os << "int11_37[" << i << "][" << j << "][" << k << "][" << l << "]=" << int11_37[i][j][k][l] << std::endl;

    // interior loops 2x2
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int m = 0; m != 4; ++m)
                for (int n = 0; n != 4; ++n)
                    for (int k = 0; k != 4; ++k)
                        for (int l = 0; l != 4; ++l)
                            os << "int22_37[" << i << "][" << j << "][" << m << "][" << l << "][" << n << "][" << k
                               << "]=" << int22_37[i][j][m][l][n][k] << std::endl;

    // interior loops 1x2
    for (int i = 0; i != 6; ++i)
        for (int j = 0; j != 6; ++j)
            for (int m = 0; m != 4; ++m)
                for (int k = 0; k != 4; ++k)
                    for (int l = 0; l != 4; ++l)
                        os << "int21_37[" << i << "][" << k << "][" << m << "][" << j << "][" << l
                           << "]=" << int21_37[i][k][m][j][l] << std::endl;

    // polyC hairpin parameters
    os << "polyC_penalty=" << polyC_penalty << std::endl
       << "polyC_slope=" << polyC_slope << std::endl
       << "polyC_int=" << polyC_int << std::endl;

    // pseudoknot energy parameters
    os << "pk_penalty=" << pk_penalty << std::endl
       << "pk_paired_penalty=" << pk_paired_penalty << std::endl
       << "pk_unpaired_penalty=" << pk_unpaired_penalty << std::endl
       << "pk_multiloop_penalty=" << pk_multiloop_penalty << std::endl
       << "pk_pk_penalty=" << pk_pk_penalty << std::endl
       << "pk_band_penalty=" << pk_band_penalty << std::endl
       << "pk_stack_span=" << pk_stack_span << std::endl
       << "pk_interior_span=" << pk_interior_span << std::endl
       << "multiloop_penalty_pk=" << multiloop_penalty_pk << std::endl
       << "multiloop_paired_penalty_pk=" << multiloop_paired_penalty_pk << std::endl
       << "multiloop_unpaired_penalty_pk" << multiloop_unpaired_penalty_pk << std::endl;

    // BIMOLECULAR TERM
    os << "intermolecular_initiation=" << intermolecular_initiation << std::endl;

    // misc
    os << "loop_greater30=" << loop_greater30 << std::endl << "hairpin_GGG=" << hairpin_GGG << std::endl;
}

float Nupack::calculate_partition_function()
{
    Q.resize(N);
    Q.fill(0.0);
    Qb.resize(N);
    Qb.fill(0.0);
    Qm.resize(N);
    Qm.fill(0.0);
    Qp.resize(N);
    Qp.fill(0.0);
    Qz.resize(N);
    Qz.fill(0.0);
    Qg.resize(N);
    Qg.fill(0.0);
    Qgl.resize(N);
    Qgl.fill(0.0);
    Qgr.resize(N);
    Qgr.fill(0.0);
    Qgls.resize(N);
    Qgls.fill(0.0);
    Qgrs.resize(N);
    Qgrs.fill(0.0);
    for (int i = 0; i != N; ++i) {
        Q(i, i - 1)  = 1.0;    // exp(- Gempty / RT) = exp(0) = 1.0
        Qz(i, i - 1) = 1.0;
    }
    DPtableX Qx, Qx1, Qx2;

    for (int l = 1; l <= N; ++l) {
        Qx.swap(Qx1);
        Qx1.swap(Qx2);
        Qx2.resize(l + 1, N);
        Qx2.fill(0.0);

        for (int i = 0; i + l <= N; ++i) {
            int j = i + l - 1;

            // Qb recursion
            if (allow_paired(i, j)) {
                Qb(i, j) = EXP(-score_hairpin(i, j) / RT);

                for (int d = i + 1; d <= j - 5; ++d)    // all possible rightmost pairs d-e
                {
                    for (int e = d + 4; e <= j - 1; ++e) {
                        if (allow_paired(d, e)) {
                            Qb(i, j) += Qb(d, e) * EXP(-score_interior(i, d, e, j, false) / RT);

                            if (d >= i + 6 && wc_pair(d, e) && wc_pair(i, j)) {
                                Qb(i, j) += Qm(i + 1, d - 1) * Qb(d, e) *
                                            EXP(
                                            -(
                                            score_multiloop(false) + score_multiloop_paired(2, false) +
                                            score_multiloop_unpaired(j - e - 1, false) + score_at_penalty(i, j) +
                                            score_at_penalty(d, e) + score_dangle(e + 1, j - 1)) /
                                            RT);
                            }
                        }
                    }
                }

                if (wc_pair(i, j)) {
                    for (int d = i + 1; d <= j - 6; ++d)    // all possible rightmost pseudoknots filling [d,e]
                    {
                        for (int e = d + 5; e <= j - 1; ++e) {
                            Qb(i, j) +=
                            Qp(d, e) * EXP(
                                       -(
                                       score_multiloop(false) + score_pk_multiloop() + score_multiloop_paired(3, false) +
                                       score_multiloop_unpaired(j - e - 1 + d - i - 1, false) + score_at_penalty(i, j) +
                                       score_dangle(e + 1, j - 1) + score_dangle(i + 1, d - 1)) /
                                       RT);

                            Qb(i, j) +=
                            Qm(i + 1, d - 1) * Qp(d, e) *
                            EXP(
                            -(
                            score_multiloop(false) + score_pk_multiloop() + score_multiloop_paired(3, false) +
                            score_multiloop_unpaired(j - e - 1, false) + score_at_penalty(i, j) + score_dangle(e + 1, j - 1)) /
                            RT);
                        }
                    }
                }
            }

            // Qg recursion
            if (allow_paired(i, j)) {
                // case 0: only 1 pair
                Qg(i, i, j, j) = 1;

                // case 1: terminal inner pair
                for (int d = i + 1; d <= j - 5; ++d) {
                    for (int e = d + 4; e <= j - 1; ++e) {
                        if (allow_paired(d, e)) Qg(i, d, e, j) += EXP(-score_interior(i, d, e, j, true) / RT);
                    }
                }
            }

            fastiloops(i, j, Qg, Qx, Qx2);

            if (allow_paired(i, j) && wc_pair(i, j)) {

                // case 2: multiloop left
                for (int d = i + 6; d <= j - 5; ++d) {
                    for (int e = d + 4; e <= j - 1; ++e) {
                        if (allow_paired(d, e) && wc_pair(d, e)) {
                            Qg(i, d, e, j) +=
                            Qm(i + 1, d - 1) *
                            EXP(
                            -(
                            score_multiloop(true) + score_multiloop_paired(2, true) + score_multiloop_unpaired(j - e - 1, true) +
                            score_at_penalty(i, j) + score_at_penalty(d, e) + score_dangle(e + 1, j - 1)) /
                            RT);
                        }
                    }
                }

                // case 3: multiloop right
                for (int d = i + 1; d <= j - 10; ++d) {
                    for (int e = d + 4; e <= j - 6; ++e) {
                        if (allow_paired(d, e) && wc_pair(d, e)) {
                            Qg(i, d, e, j) +=
                            Qm(e + 1, j - 1) *
                            EXP(
                            -(
                            score_multiloop(true) + score_multiloop_paired(2, true) + score_multiloop_unpaired(d - i - 1, true) +
                            score_at_penalty(i, j) + score_at_penalty(d, e) + score_dangle(i + 1, d - 1)) /
                            RT);
                        }
                    }
                }

                // case 4: multiloop both sides
                for (int d = i + 6; d <= j - 10; ++d) {
                    for (int e = d + 4; e <= j - 6; ++e) {
                        if (allow_paired(d, e) && wc_pair(d, e)) {
                            Qg(i, d, e, j) +=
                            Qm(i + 1, d - 1) * Qm(e + 1, j - 1) *
                            EXP(-(score_multiloop(true) + score_multiloop_paired(2, true) + score_at_penalty(i, j) + score_at_penalty(d, e)) / RT);
                        }
                    }
                }

                // case 5: interior loop + multi left
                for (int d = i + 7; d <= j - 6; ++d) {
                    for (int e = d + 4; e <= j - 2; ++e) {
                        if (allow_paired(d, e)) {
                            for (int f = e + 1; f <= j - 1; ++f) {
                                Qg(i, d, e, j) +=
                                Qgls(i + 1, d, e, f) * EXP(
                                                       -(
                                                       score_multiloop(true) + score_multiloop_paired(1, true) +
                                                       score_multiloop_unpaired(j - f - 1, true) +
                                                       score_at_penalty(i, j) + score_dangle(f + 1, j - 1)) /
                                                       RT);
                            }
                        }
                    }
                }

                // case 6: interior loop + multi right
                for (int d = i + 2; d <= j - 11; ++d) {
                    for (int e = d + 4; e <= j - 7; ++e) {
                        if (allow_paired(d, e)) {
                            for (int c = i + 1; c <= d - 1; ++c) {
                                Qg(i, d, e, j) +=
                                Qgrs(c, d, e, j - 1) * EXP(
                                                       -(
                                                       score_multiloop(true) + score_multiloop_paired(1, true) +
                                                       score_multiloop_unpaired(c - i - 1, true) +
                                                       score_at_penalty(i, j) + score_dangle(i + 1, c - 1)) /
                                                       RT);
                            }
                        }
                    }
                }

                // case 7: interior loop + multi both sides
                for (int d = i + 7; d <= j - 11; ++d) {
                    for (int e = d + 4; e <= j - 7; ++e) {
                        if (allow_paired(d, e)) {
                            for (int c = i + 6; c <= d - 1; ++c) {
                                Qg(i, d, e, j) +=
                                Qm(i + 1, c - 1) * Qgrs(c, d, e, j - 1) *
                                EXP(-(score_multiloop(true) + score_multiloop_paired(1, true) + score_at_penalty(i, j)) / RT);
                            }
                        }
                    }
                }
            }

            // Qgls recursion
            for (int c = i + 5; c <= j - 6; ++c) {
                if (allow_paired(c, j) && wc_pair(c, j)) {
                    for (int d = c + 1; d <= j - 5; ++d) {
                        for (int e = d + 4; e <= j - 1; ++e) {
                            if (allow_paired(d, e)) {
                                Qgls(i, d, e, j) += Qm(i, c - 1) * Qg(c, d, e, j) *
                                                    EXP(-(score_multiloop_paired(1, true) + score_at_penalty(c, j)) / RT);
                            }
                        }
                    }
                }
            }

            // Qgrs recursion
            for (int d = i + 1; d <= j - 10; ++d) {
                for (int e = d + 4; e <= j - 6; ++e) {
                    if (allow_paired(d, e)) {
                        for (int f = e + 1; f <= j - 5; ++f) {
                            if (allow_paired(i, f) && wc_pair(i, f)) {
                                Qgrs(i, d, e, j) += Qg(i, d, e, f) * Qm(f + 1, j) *
                                                    EXP(-(score_multiloop_paired(1, true) + score_at_penalty(i, f)) / RT);
                            }
                        }
                    }
                }
            }

            // Qgl recursions
            for (int d = i + 1; d <= j - 5; ++d) {
                for (int f = d + 4; f <= j - 1; ++f) {
                    if (allow_paired(d, f) && wc_pair(d, f)) {
                        for (int e = d; e <= f - 2; ++e)    // f-3???
                        {
                            Qgl(i, e, f, j) +=
                            Qg(i, d, f, j) * Qz(d + 1, e) * EXP(-(score_pk_paired(1) + score_at_penalty(d, f)) / RT);
                        }
                    }
                }
            }

            // Qgr recursion
            for (int d = i + 1; d <= j - 3; ++d) {
                for (int e = d + 2; e <= j - 1; ++e) {
                    for (int f = e; f <= j - 1; ++f) {
                        Qgr(i, d, e, j) += Qgl(i, d, f, j) * Qz(e, f - 1);
                    }
                }
            }

            // Qp recursion
            // case 1: both Qg are exactly 1 pair
            // first case is exactly 1 pair per Og
            if (j - i > 4) {
                int a = i;
                int f = j;
                for (int b = a + 1; b <= j - 4; ++b) {
                    if (allow_paired(b, j) && wc_pair(b, j)) {
                        int c = b;
                        for (int d = std::max(c + 1, a + 4); d <= j - 1; ++d) {
                            if (allow_paired(a, d) && wc_pair(a, d)) {
                                int e = d;
                                Qp(i, j) += Qg(i, a, d, e) * Qg(b, c, f, j) * Qz(e + 1, f - 1) * Qz(c + 1, d - 1) *
                                            Qz(a + 1, b - 1) *
                                            EXP(
                                            -(
                                            score_pk_paired(2) + score_pk_band(2) + score_at_penalty(a, d) +
                                            score_at_penalty(c, f) + score_at_penalty(i, e) + score_at_penalty(b, j)) /
                                            RT);
                            }
                        }
                    }
                }
            }

            if (j - i > 6) {
                // case 2 left Og is exactly 1 pair, right is 2+
                for (int d = i + 1; d <= j - 6; ++d) {
                    if (allow_paired(d, j) && wc_pair(d, j)) {
                        for (int e = std::max(d + 2, i + 4); e <= j - 2; ++e) {
                            int f = e;
                            if (allow_paired(i, f) && wc_pair(i, f)) {
                                Qp(i, j) +=
                                Qg(i, i, e, f) * Qz(i + 1, d - 1) * Qgr(d, e - 1, f + 1, j) *
                                EXP(-(score_pk_paired(1) + score_pk_band(2) + score_at_penalty(d, j) + score_at_penalty(i, f) * 2) / RT);
                            }
                        }
                    }
                }

                // case 2 left Qg is 2+ pairs, right is 1
                for (int d = i + 2; d <= j - 4; ++d) {
                    if (allow_paired(d, j) && wc_pair(d, j)) {
                        for (int e = std::max(d + 1, i + 4); e <= j - 2; ++e) {
                            for (int f = e + 1; f <= j - 1; ++f) {
                                if (allow_paired(i, f) && wc_pair(i, f)) {
                                    Qp(i, j) +=
                                    Qgl(i, d - 1, e, f) * Qg(d, d, j, j) * Qz(d + 1, e - 1) * Qz(f + 1, j - 1) *
                                    EXP(-(score_pk_paired(1) + score_pk_band(2) + score_at_penalty(d, j) * 2 + score_at_penalty(i, f)) / RT);
                                }
                            }
                        }
                    }
                }
            }

            // otherwise
            if (j - i > 7) {
                for (int d = i + 2; d <= j - 4; ++d) {
                    if (allow_paired(d, j) && wc_pair(d, j)) {
                        for (int e = std::max(d + 2, i + 5); e <= j - 3; ++e) {
                            for (int f = e + 1; f <= j - 2; ++f) {
                                if (allow_paired(i, f) && wc_pair(i, f)) {
                                    Qp(i, j) += Qgl(i, d - 1, e, f) * Qgr(d, e - 1, f + 1, j) *
                                                EXP(-(score_pk_band(2) + score_at_penalty(d, j) + score_at_penalty(i, j)) / RT);
                                }
                            }
                        }
                    }
                }
            }

            // Q, Qm, Qz recursions
            Q(i, j) = EXP(-score_dangle(i, j) / RT);    // empty recursion

            if (i != 0 && j != N - 1) {
                Qz(i, j) = EXP(-(score_dangle(i, j) + score_pk_unpaired(j - i + 1)) / RT);
            }

            for (int d = i; d <= j - 4; ++d) {
                for (int e = d + 4; e <= j; ++e) {
                    if (allow_paired(d, e) && wc_pair(d, e)) {
                        Q(i, j) += Q(i, d - 1) * Qb(d, e) * EXP(-(score_at_penalty(d, e) + score_dangle(e + 1, j)) / RT);

                        if (i != 0 && j != N - 1) {
                            Qm(i, j) +=
                            Qb(d, e) * EXP(
                                       -(
                                       score_multiloop_paired(1, false) + score_multiloop_unpaired(d - i + j - e, false) +
                                       score_at_penalty(d, e) + score_dangle(e + 1, j) + score_dangle(i, d - 1)) /
                                       RT);

                            if (d >= i + 5) {
                                Qm(i, j) += Qm(i, d - 1) * Qb(d, e) *
                                            EXP(
                                            -(
                                            score_multiloop_paired(1, false) + score_multiloop_unpaired(j - e, false) +
                                            score_at_penalty(d, e) + score_dangle(e + 1, j)) /
                                            RT);
                            }

                            Qz(i, j) +=
                            Qz(i, d - 1) * Qb(d, e) *
                            EXP(-(score_pk_paired(1) + score_pk_unpaired(j - e) + score_at_penalty(d, e) + score_dangle(e + 1, j)) / RT);
                        }
                    }
                }
            }

            for (int d = i; d <= j - 5; ++d) {
                for (int e = d + 5; e <= j; ++e) {
                    Q(i, j) += Q(i, d - 1) * Qp(d, e) * EXP(-(score_pk() + score_dangle(e + 1, j)) / RT);

                    if (i != 0 && j != N - 1) {
                        Qm(i, j) += Qp(d, e) * EXP(
                                               -(
                                               score_pk_multiloop() + score_multiloop_paired(2, false) +
                                               score_multiloop_unpaired(d - i + j - e, false) + score_dangle(e + 1, j) +
                                               score_dangle(i, d - 1)) /
                                               RT);

                        if (d >= i + 5) {
                            Qm(i, j) += Qm(i, d - 1) * Qp(d, e) *
                                        EXP(
                                        -(
                                        score_pk_multiloop() + score_multiloop_paired(2, false) +
                                        score_multiloop_unpaired(j - e, false) + score_dangle(e + 1, j)) /
                                        RT);
                        }

                        Qz(i, j) +=
                        Qz(i, d - 1) * Qp(d, e) *
                        EXP(-(score_pk_pk() + score_pk_paired(2) + score_pk_unpaired(j - e) + score_dangle(e + 1, j)) / RT);
                    }
                }
            }
            printf("%d,%d: %Lf\n", i, j, Q(i, j));
        }
    }

    return Q(0, N - 1);
}

void Nupack::fastiloops(int i, int j, DPtable4& Qg, DPtableX& Qx, DPtableX& Qx2)
{
    int l = j - i + 1;
    if (l >= 17)    // smallest subsequence not added to Qg as special case
    {
        for (int d = i + 6; d <= j - 10; ++d) {
            for (int e = d + 4; e <= j - 6; ++e) {
                if (allow_paired(d, e)) {
                    int l1 = 4;    // explicitly add in terms for l1=4, l2>=4
                    int c  = i + l1 + 1;
                    for (int l2 = 4; l2 <= j - e - 2; ++l2) {
                        int s = l1 + l2;
                        int f = j - l2 - 1;
                        if (allow_paired(c, f)) {
                            Qx(i, d, e, s) +=
                            Qg(c, d, e, f) *
                            EXP(-(score_interior_asymmetry(l1, l2) + score_interior_mismatch(f, c, f + 1, c - 1)) / RT);
                        }
                    }

                    if (d >= i + 7) {
                        int l2 = 4;    // explicitly add in terms of l1>=5, l2=4
                        int f  = j - l2 - 1;
                        for (int l1 = 5; l1 <= d - i - 2; ++l1) {
                            int s = l1 + l2;
                            int c = i + l1 + 1;
                            if (allow_paired(c, f)) {
                                Qx(i, d, e, s) +=
                                Qg(c, d, e, f) *
                                EXP(-(score_interior_asymmetry(l1, l2) + score_interior_mismatch(f, c, f + 1, c - 1)) / RT);
                            }
                        }
                    }
                }
            }
        }
    }

    for (int d = i + 1; d <= j - 5; ++d) {
        for (int e = d + 4; e <= j - 1; ++e) {
            if (allow_paired(d, e)) {
                // conveRT Qx into interior loop energies
                if (l >= 17 && allow_paired(i, j)) {
                    for (int s = 8; s <= l - 9; ++s) {
                        Qg(i, d, e, j) += Qx(i, d, e, s) * EXP(-score_interior_mismatch(i, j, i + 1, j - 1) / RT);
                    }
                }

                // extend loops for future use
                if (i != 0 && j != N - 1) {
                    for (int s = 8; s <= l - 9; ++s) {
                        Qx2(i - 1, d, e, s + 2) = Qx(i, d, e, s) * EXP(-(score_loop(s + 2) - score_loop(s)) / RT);
                    }
                }

                if (allow_paired(i, j)) {
                    // Add small inextensible interior loops to Qg as special cases
                    for (int l1 = 0; l1 <= std::min(3, d - i - 2); ++l1) {
                        int c = i + l1 + 1;
                        for (int l2 = 0; l2 <= std::min(3, j - e - 2); ++l2) {
                            int f = j - l2 - 1;
                            if (allow_paired(c, f)) {
                                Qg(i, d, e, j) += Qg(c, d, e, f) * EXP(-score_interior(i, c, f, j, true) / RT);
                            }
                        }
                    }
                    // Add bulge loops and large asymmetric loops as special cases
                    for (int l1 = 0; l1 <= std::min(3, d - i - 2); ++l1)    // cases l1=0,1,2,3, l2>=4
                    {
                        int c = i + l1 + 1;
                        for (int l2 = 4; l2 <= j - e - 2; ++l2) {
                            int f = j - l2 - 1;
                            if (allow_paired(c, f)) {
                                Qg(i, d, e, j) += Qg(c, d, e, f) * EXP(-score_interior(i, c, f, j, true) / RT);
                            }
                        }
                    }
                    for (int l2 = 0; l2 <= std::min(3, j - e - 2); ++l2) {
                        int f = j - l2 - 1;
                        for (int l1 = 4; l1 <= d - i - 2; ++l1) {
                            int c = i + l1 + 1;
                            if (allow_paired(c, f)) {
                                Qg(i, d, e, j) += Qg(c, d, e, f) * EXP(-score_interior(i, c, f, j, true) / RT);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Nupack::calculate_posterior()
{
    // std::cout << "calculate_posterior_nupack\n";
    P.resize(N);
    P.fill(0.0);
    Pb.resize(N);
    Pb.fill(0.0);
    Pm.resize(N);
    Pm.fill(0.0);
    Pp.resize(N);
    Pp.fill(0.0);
    Pz.resize(N);
    Pz.fill(0.0);
    Pbg.resize(N);
    Pbg.fill(0.0);
    Pg.resize(N);
    Pg.fill(0.0);
    Pgl.resize(N);
    Pgl.fill(0.0);
    Pgr.resize(N);
    Pgr.fill(0.0);
    Pgls.resize(N);
    Pgls.fill(0.0);
    Pgrs.resize(N);
    Pgrs.fill(0.0);
    P(0, N - 1) = 1.0;
    DPtableX Qx, Qx1, Qx2;
    Qx.resize(N, N);
    Qx1.resize(N - 1, N);
    Qx2.resize(N - 2, N);
    DPtableX Px, Px1, Px2;
    Px.resize(N, N);
    Px1.resize(N - 1, N);
    Px2.resize(N - 2, N);

    for (int l = N; l >= 1; --l) {
        Qx.swap(Qx1);
        Qx1.swap(Qx2);
        Qx2.resize(l - 3, N);
        Qx2.fill(0.0);
        Px.swap(Px1);
        Px1.swap(Px2);
        Px2.resize(l - 3, N);
        Px2.fill(0.0);

        for (int i = 0; i + l <= N; ++i) {
            int j = i + l - 1;

            // P, Pm, Pz recursions
            for (int d = i; d <= j - 4; ++d) {
                for (int e = d + 4; e <= j; ++e) {
                    if (allow_paired(d, e) && wc_pair(d, e)) {
                        if (P(i, j) > 0.0) {
                            float p = Q(i, d - 1) * Qb(d, e) *
                                      EXP(-(score_at_penalty(d, e) + score_dangle(e + 1, j)) / RT) / Q(i, j) * P(i, j);
                            P(i, d - 1) += p;
                            Pb(d, e) += p;
                            assert(!std::isnan(p));
                        }

                        if (i != 0 && j != N - 1) {
                            if (Pm(i, j) > 0.0) {
                                float p = Qb(d, e) *
                                          EXP(
                                          -(
                                          score_multiloop_paired(1, false) + score_multiloop_unpaired(d - i + j - e, false) +
                                          score_at_penalty(d, e) + score_dangle(e + 1, j) + score_dangle(i, d - 1)) /
                                          RT) /
                                          Qm(i, j) * Pm(i, j);
                                Pb(d, e) += p;

                                if (d >= i + 5) {
                                    p = Qm(i, d - 1) * Qb(d, e) *
                                        EXP(
                                        -(
                                        score_multiloop_paired(1, false) + score_multiloop_unpaired(j - e, false) +
                                        score_at_penalty(d, e) + score_dangle(e + 1, j)) /
                                        RT) /
                                        Qm(i, j) * Pm(i, j);
                                    Pm(i, d - 1) += p;
                                    Pb(d, e) += p;
                                    assert(!std::isnan(p));
                                }
                            }

                            if (Pz(i, j) > 0.0) {
                                float p =
                                Qz(i, d - 1) * Qb(d, e) *
                                EXP(-(score_pk_paired(1) + score_pk_unpaired(j - e) + score_at_penalty(d, e) + score_dangle(e + 1, j)) / RT) /
                                Qz(i, j) * Pz(i, j);
                                Pz(i, d - 1) += p;
                                Pb(d, e) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }

            for (int d = i; d <= j - 5; ++d) {
                for (int e = d + 5; e <= j; ++e) {
                    if (P(i, j) > 0.0) {
                        float p = Q(i, d - 1) * Qp(d, e) * EXP(-(score_pk() + score_dangle(e + 1, j)) / RT) / Q(i, j) * P(i, j);
                        P(i, d - 1) += p;
                        Pp(d, e) += p;
                        assert(!std::isnan(p));
                    }

                    if (i != 0 && j != N - 1) {
                        if (Pm(i, j) > 0.0) {
                            float p =
                            Qp(d, e) *
                            EXP(
                            -(
                            score_pk_multiloop() + score_multiloop_paired(2, false) +
                            score_multiloop_unpaired(d - i + j - e, false) + score_dangle(e + 1, j) + score_dangle(i, d - 1)) /
                            RT) /
                            Qm(i, j) * Pm(i, j);
                            Pp(d, e) += p;
                            assert(!std::isnan(p));

                            if (d >= i + 5) {
                                p = Qm(i, d - 1) * Qp(d, e) *
                                    EXP(
                                    -(
                                    score_pk_multiloop() + score_multiloop_paired(2, false) +
                                    score_multiloop_unpaired(j - e, false) + score_dangle(e + 1, j)) /
                                    RT) /
                                    Qm(i, j) * Pm(i, j);
                                Pm(i, d - 1) += p;
                                Pp(d, e) += p;
                                assert(!std::isnan(p));
                            }
                        }

                        if (Pz(i, j) > 0.0) {
                            float p =
                            Qz(i, d - 1) * Qp(d, e) *
                            EXP(-(score_pk_pk() + score_pk_paired(2) + score_pk_unpaired(j - e) + score_dangle(e + 1, j)) / RT) /
                            Qz(i, j) * Pz(i, j);
                            Pz(i, d - 1) += p;
                            Pp(d, e) += p;
                            assert(!std::isnan(p));
                        }
                    }
                }
            }

            // Pp recursion
            // case 1: both Qg are exactly 1 pair
            // first case is exactly 1 pair per Og
            if (j - i > 4) {
                int a = i;
                int f = j;
                for (int b = a + 1; b <= j - 4; ++b) {
                    if (allow_paired(b, j) && wc_pair(b, j)) {
                        int c = b;
                        for (int d = std::max(c + 1, a + 4); d <= j - 1; ++d) {
                            if (allow_paired(a, d) && wc_pair(a, d) && Pp(i, j) > 0.0) {
                                int   e = d;
                                float p = Qg(i, a, d, e) * Qg(b, c, f, j) * Qz(e + 1, f - 1) * Qz(c + 1, d - 1) *
                                          Qz(a + 1, b - 1) *
                                          EXP(
                                          -(
                                          score_pk_paired(2) + score_at_penalty(a, d) + score_at_penalty(c, f) +
                                          score_at_penalty(i, e) + score_at_penalty(b, j)) /
                                          RT) /
                                          Qp(i, j) * Pp(i, j);
                                Pg(i, a, d, e) += p;
                                Pg(b, c, f, j) += p;
                                Pz(e + 1, f - 1) += p;
                                Pz(c + 1, d - 1) += p;
                                Pz(a + 1, b - 1) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }

            if (j - i > 6) {
                // case 2 left Og is exactly 1 pair, right is 2+
                for (int d = i + 1; d <= j - 6; ++d) {
                    if (allow_paired(d, j) && wc_pair(d, j)) {
                        for (int e = std::max(d + 2, i + 4); e <= j - 2; ++e) {
                            int f = e;
                            if (allow_paired(i, f) && wc_pair(i, f) && Pp(i, j) > 0.0) {
                                float p = Qg(i, i, e, f) * Qz(i + 1, d - 1) * Qgr(d, e - 1, f + 1, j) *
                                          EXP(-(score_pk_paired(1) + score_at_penalty(d, j) + score_at_penalty(i, f) * 2) / RT) /
                                          Qp(i, j) * Pp(i, j);
                                Pg(i, i, e, f) += p;
                                Pz(i + 1, d - 1) += p;
                                Pgr(d, e - 1, f + 1, j) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }

                // case 2 left Qg is 2+ pairs, right is 1
                for (int d = i + 2; d <= j - 4; ++d) {
                    if (allow_paired(d, j) && wc_pair(d, j)) {
                        for (int e = std::max(d + 1, i + 4); e <= j - 2; ++e) {
                            for (int f = e + 1; f <= j - 1; ++f) {
                                if (allow_paired(i, f) && wc_pair(i, f) && Pp(i, j) > 0.0) {
                                    float p =
                                    Qgl(i, d - 1, e, f) * Qg(d, d, j, j) * Qz(d + 1, e - 1) * Qz(f + 1, j - 1) *
                                    EXP(-(score_pk_paired(1) + score_at_penalty(d, j) * 2 + score_at_penalty(i, f)) / RT) /
                                    Qp(i, j) * Pp(i, j);
                                    Pgl(i, d - 1, e, f) += p;
                                    Pg(d, d, j, j) += p;
                                    Pz(d + 1, e - 1) += p;
                                    Pz(f + 1, j - 1) += p;
                                    assert(!std::isnan(p));
                                }
                            }
                        }
                    }
                }
            }

            // otherwise
            if (j - i > 7) {
                for (int d = i + 2; d <= j - 4; ++d) {
                    if (allow_paired(d, j) && wc_pair(d, j)) {
                        for (int e = std::max(d + 2, i + 5); e <= j - 3; ++e) {
                            for (int f = e + 1; f <= j - 2; ++f) {
                                if (allow_paired(i, f) && wc_pair(i, f) && Pp(i, j) > 0.0) {
                                    float p = Qgl(i, d - 1, e, f) * Qgr(d, e - 1, f + 1, j) *
                                              EXP(-(score_at_penalty(d, j) + score_at_penalty(i, j)) / RT) / Qp(i, j) *
                                              Pp(i, j);
                                    Pgl(i, d - 1, e, f) += p;
                                    Pgr(d, e - 1, f + 1, j) += p;
                                    assert(!std::isnan(p));
                                }
                            }
                        }
                    }
                }
            }

            // Pgr recursion
            for (int d = i + 1; d <= j - 3; ++d) {
                for (int e = d + 2; e <= j - 1; ++e) {
                    for (int f = e; f <= j - 1; ++f) {
                        if (Pgr(i, d, e, j) > 0.0) {
                            float p = Qgl(i, d, f, j) * Qz(e, f - 1) / Qgr(i, d, e, j) * Pgr(i, d, e, j);
                            Pgl(i, d, f, j) += p;
                            Pz(e, f - 1) += p;
                            assert(!std::isnan(p));
                        }
                    }
                }
            }

            // Pgl recursions
            for (int d = i + 1; d <= j - 5; ++d) {
                for (int f = d + 4; f <= j - 1; ++f) {
                    if (allow_paired(d, f) && wc_pair(d, f)) {
                        for (int e = d; e <= f - 2; ++e)    // f-3???
                        {
                            if (Pgl(i, e, f, j) > 0.0) {
                                float p = Qg(i, d, f, j) * Qz(d + 1, e) *
                                          EXP(-(score_pk_paired(1) + score_at_penalty(d, f)) / RT) / Qgl(i, e, f, j) *
                                          Pgl(i, e, f, j);
                                Pg(i, d, f, j) += p;
                                Pz(d + 1, e) += p;
                                Pbg(d, f) += p;    // Pbg inner gap-spanning base-pairing prob
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }

            // Pgrs recursion
            for (int d = i + 1; d <= j - 10; ++d) {
                for (int e = d + 4; e <= j - 6; ++e) {
                    if (allow_paired(d, e)) {
                        for (int f = e + 1; f <= j - 5; ++f) {
                            if (allow_paired(i, f) && wc_pair(i, f) && Pgrs(i, d, e, j) > 0.0) {
                                float p = Qg(i, d, e, f) * Qm(f + 1, j) *
                                          EXP(-(score_multiloop_paired(1, true) + score_at_penalty(i, f)) / RT) /
                                          Qgrs(i, d, e, j) * Pgrs(i, d, e, j);
                                Pg(i, d, e, f) += p;
                                Pm(f + 1, j) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }

            // Pgls recursion
            for (int c = i + 5; c <= j - 6; ++c) {
                if (allow_paired(c, j) && wc_pair(c, j)) {
                    for (int d = c + 1; d <= j - 5; ++d) {
                        for (int e = d + 4; e <= j - 1; ++e) {
                            if (allow_paired(d, e) && Pgls(i, d, e, j) > 0.0) {
                                float p = Qm(i, c - 1) * Qg(c, d, e, j) *
                                          EXP(-(score_multiloop_paired(1, true) + score_at_penalty(c, j)) / RT) /
                                          Qgls(i, d, e, j) * Pgls(i, d, e, j);
                                Pm(i, c - 1) += p;
                                Pg(c, d, e, j) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }

            // Pg recursion
            fastiloops_pr(i, j, Qg, Qx, Qx2, Pg, Px, Px2);

            if (allow_paired(i, j) && wc_pair(i, j)) {
                // case 2: multiloop left
                for (int d = i + 6; d <= j - 5; ++d) {
                    for (int e = d + 4; e <= j - 1; ++e) {
                        if (allow_paired(d, e) && wc_pair(d, e) && Pg(i, d, e, j) > 0.0) {
                            float p =
                            Qm(i + 1, d - 1) *
                            EXP(
                            -(
                            score_multiloop(true) + score_multiloop_paired(2, true) + score_multiloop_unpaired(j - e - 1, true) +
                            score_at_penalty(i, j) + score_at_penalty(d, e) + score_dangle(e + 1, j - 1)) /
                            RT) /
                            Qg(i, d, e, j) * Pg(i, d, e, j);
                            Pm(i + 1, d - 1) += p;
                            assert(!std::isnan(p));
                        }
                    }
                }

                // case 3: multiloop right
                for (int d = i + 1; d <= j - 10; ++d) {
                    for (int e = d + 4; e <= j - 6; ++e) {
                        if (allow_paired(d, e) && wc_pair(d, e) && Pg(i, d, e, j) > 0.0) {
                            float p =
                            Qm(e + 1, j - 1) *
                            EXP(
                            -(
                            score_multiloop(true) + score_multiloop_paired(2, true) + score_multiloop_unpaired(d - i - 1, true) +
                            score_at_penalty(i, j) + score_at_penalty(d, e) + score_dangle(i + 1, d - 1)) /
                            RT) /
                            Qg(i, d, e, j) * Pg(i, d, e, j);
                            Pm(e + 1, j - 1) += p;
                            assert(!std::isnan(p));
                        }
                    }
                }

                // case 4: multiloop both sides
                for (int d = i + 6; d <= j - 10; ++d) {
                    for (int e = d + 4; e <= j - 6; ++e) {
                        if (allow_paired(d, e) && wc_pair(d, e) && Pg(i, d, e, j) > 0.0) {
                            float p =
                            Qm(i + 1, d - 1) * Qm(e + 1, j - 1) *
                            EXP(-(score_multiloop(true) + score_multiloop_paired(2, true) + score_at_penalty(i, j) + score_at_penalty(d, e)) / RT) /
                            Qg(i, d, e, j) * Pg(i, d, e, j);
                            Pm(i + 1, d - 1) += p;
                            Pm(e + 1, j - 1) += p;
                            assert(!std::isnan(p));
                        }
                    }
                }

                // case 5: interior loop + multi left
                for (int d = i + 7; d <= j - 6; ++d) {
                    for (int e = d + 4; e <= j - 2; ++e) {
                        if (allow_paired(d, e)) {
                            for (int f = e + 1; f <= j - 1; ++f) {
                                if (Pg(i, d, e, j) > 0.0) {
                                    float p = Qgls(i + 1, d, e, f) *
                                              EXP(
                                              -(
                                              score_multiloop(true) + score_multiloop_paired(1, true) +
                                              score_multiloop_unpaired(j - f - 1, true) + score_at_penalty(i, j) +
                                              score_dangle(f + 1, j - 1)) /
                                              RT) /
                                              Qg(i, d, e, j) * Pg(i, d, e, j);
                                    Pgls(i + 1, d, e, f) += p;
                                    assert(!std::isnan(p));
                                }
                            }
                        }
                    }
                }

                // case 6: interior loop + multi right
                for (int d = i + 2; d <= j - 11; ++d) {
                    for (int e = d + 4; e <= j - 7; ++e) {
                        if (allow_paired(d, e)) {
                            for (int c = i + 1; c <= d - 1; ++c) {
                                if (Pg(i, d, e, j) > 0.0) {
                                    float p = Qgrs(c, d, e, j - 1) *
                                              EXP(
                                              -(
                                              score_multiloop(true) + score_multiloop_paired(1, true) +
                                              score_multiloop_unpaired(c - i - 1, true) + score_at_penalty(i, j) +
                                              score_dangle(i + 1, c - 1)) /
                                              RT) /
                                              Qg(i, d, e, j) * Pg(i, d, e, j);
                                    Pgrs(c, d, e, j - 1) += p;
                                    assert(!std::isnan(p));
                                }
                            }
                        }
                    }
                }

                // case 7: interior loop + multi both sides
                for (int d = i + 7; d <= j - 11; ++d) {
                    for (int e = d + 4; e <= j - 7; ++e) {
                        if (allow_paired(d, e)) {
                            for (int c = i + 6; c <= d - 1; ++c) {
                                if (Pg(i, d, e, j) > 0.0) {
                                    float p =
                                    Qm(i + 1, c - 1) * Qgrs(c, d, e, j - 1) *
                                    EXP(-(score_multiloop(true) + score_multiloop_paired(1, true) + score_at_penalty(i, j)) / RT) /
                                    Qg(i, d, e, j) * Pg(i, d, e, j);
                                    Pm(i + 1, c - 1) += p;
                                    Pgrs(c, d, e, j - 1) += p;
                                    assert(!std::isnan(p));
                                }
                            }
                        }
                    }
                }
            }

            // Pbg outer gap-spanning base-pairing prob
            for (int d = i + 1; d <= j - 5; ++d)
                for (int e = d + 4; e <= j - 1; ++e) Pbg(i, j) += Pg(i, d, e, j);

            // Pb recursion
            if (allow_paired(i, j)) {
                for (int d = i + 1; d <= j - 5; ++d)    // all possible rightmost pairs d-e
                {
                    for (int e = d + 4; e <= j - 1; ++e) {
                        if (allow_paired(d, e)) {
                            if (Pb(i, j) > 0.0) {
                                float p = Qb(d, e) * EXP(-score_interior(i, d, e, j, false) / RT) / Qb(i, j) * Pb(i, j);
                                Pb(d, e) += p;
                                assert(!std::isnan(p));

                                if (d >= i + 6 && wc_pair(d, e) && wc_pair(i, j)) {
                                    float p = Qm(i + 1, d - 1) * Qb(d, e) *
                                              EXP(
                                              -(
                                              score_multiloop(false) + score_multiloop_paired(2, false) +
                                              score_multiloop_unpaired(j - e - 1, false) + score_at_penalty(i, j) +
                                              score_at_penalty(d, e) + score_dangle(e + 1, j - 1)) /
                                              RT) /
                                              Qb(i, j) * Pb(i, j);
                                    Pm(i + 1, d - 1) += p;
                                    Pb(d, e) += p;
                                    assert(!std::isnan(p));
                                }
                            }
                        }
                    }
                }

                if (wc_pair(i, j)) {
                    for (int d = i + 1; d <= j - 6; ++d)    // all possible rightmost pseudoknots filling [d,e]
                    {
                        for (int e = d + 5; e <= j - 1; ++e) {
                            if (Pb(i, j) > 0.0) {
                                float p;
                                p = Qp(d, e) *
                                    EXP(
                                    -(
                                    score_multiloop(false) + score_pk_multiloop() + score_multiloop_paired(3, false) +
                                    score_multiloop_unpaired(j - e - 1 + d - i - 1, false) + score_at_penalty(i, j) +
                                    score_dangle(e + 1, j - 1) + score_dangle(i + 1, d - 1)) /
                                    RT) /
                                    Qb(i, j) * Pb(i, j);
                                Pp(d, e) += p;
                                assert(!std::isnan(p));

                                p = Qm(i + 1, d - 1) * Qp(d, e) *
                                    EXP(
                                    -(
                                    score_multiloop(false) + score_pk_multiloop() + score_multiloop_paired(3, false) +
                                    score_multiloop_unpaired(j - e - 1, false) + score_at_penalty(i, j) +
                                    score_dangle(e + 1, j - 1)) /
                                    RT) /
                                    Qb(i, j) * Pb(i, j);
                                Pm(i + 1, d - 1) += p;
                                Pp(d, e) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }
        }
    }
}

void Nupack::fastiloops_pr(int i, int j, DPtable4& Qg, DPtableX& Qx, DPtableX& Qx2, DPtable4& Pg, DPtableX& Px, DPtableX& Px2)
{
    int l = j - i + 1;

    if (allow_paired(i, j)) {
        for (int d = i + 1; d <= j - 5; ++d) {
            for (int e = d + 4; e <= j - 1; ++e) {
                if (allow_paired(d, e)) {
                    // Add small inextensible interior loops to Qg as special cases
                    for (int l1 = 0; l1 <= std::min(3, d - i - 2); ++l1) {
                        int c = i + l1 + 1;
                        for (int l2 = 0; l2 <= std::min(3, j - e - 2); ++l2) {
                            int f = j - l2 - 1;
                            if (allow_paired(c, f) && Pg(i, d, e, j) > 0.0) {
                                float p = Qg(c, d, e, f) * EXP(-score_interior(i, c, f, j, true) / RT) /
                                          Qg(i, d, e, j) * Pg(i, d, e, j);
                                Pg(c, d, e, f) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                    // Add bulge loops and large asymmetric loops as special cases
                    for (int l1 = 0; l1 <= std::min(3, d - i - 2); ++l1)    // cases l1=0,1,2,3, l2>=4
                    {
                        int c = i + l1 + 1;
                        for (int l2 = 4; l2 <= j - e - 2; ++l2) {
                            int f = j - l2 - 1;
                            if (allow_paired(c, f) && Pg(i, d, e, j) > 0.0) {
                                float p = Qg(c, d, e, f) * EXP(-score_interior(i, c, f, j, true) / RT) /
                                          Qg(i, d, e, j) * Pg(i, d, e, j);
                                Pg(c, d, e, f) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                    for (int l2 = 0; l2 <= std::min(3, j - e - 2); ++l2) {
                        int f = j - l2 - 1;
                        for (int l1 = 4; l1 <= d - i - 2; ++l1) {
                            int c = i + l1 + 1;
                            if (allow_paired(c, f) && Pg(i, d, e, j) > 0.0) {
                                float p = Qg(c, d, e, f) * EXP(-score_interior(i, c, f, j, true) / RT) /
                                          Qg(i, d, e, j) * Pg(i, d, e, j);
                                Pg(c, d, e, f) += p;
                                assert(!std::isnan(p));
                            }
                        }
                    }
                }
            }
        }
    }

    // Add cases that are at an end with l1>=4, l2>=4
    if ((i == 0 || j == N - 1) && l >= 17) {
        for (int d = i + 6; d <= j - 10; ++d) {
            for (int e = d + 4; e <= j - 6; ++e) {
                if (allow_paired(d, e)) {
                    for (int c = i + 5; c <= d - 1; ++c) {
                        for (int f = e + 1; f <= j - 5; ++f) {
                            if (allow_paired(c, f)) {
                                int l1 = c - i - 1;
                                int l2 = j - f - 1;
                                int s  = l1 + l2;
                                Qx(i, d, e, s) +=
                                Qg(c, d, e, f) *
                                EXP(-(score_interior_asymmetry(l1, l2) + score_interior_mismatch(f, c, f + 1, c - 1)) / RT);
                            }
                        }
                    }
                }
            }
        }
    }

    // Use Qx to finish calculation of Px
    if (allow_paired(i, j)) {
        for (int d = i + 1; d <= j - 5; ++d) {
            for (int e = d + 4; e <= j - 1; ++e) {
                if (allow_paired(d, e)) {
                    for (int s = 8; s <= l - 9; ++s) {
                        float p = Qx(i, d, e, s) * EXP(-score_interior_mismatch(i, j, i + 1, j - 1) / RT) /
                                  Qg(i, d, e, j) * Pg(i, d, e, j);
                        Px(i, d, e, s) += p;
                        assert(!std::isnan(p));
                    }
                }
            }
        }
    }

    // Calculate Pg contribution using Qx and Px
    if (l >= 17) {
        for (int d = i + 6; d <= j - 10; ++d) {
            for (int e = d + 4; e <= j - 6; ++e) {
                if (allow_paired(d, e)) {
                    int l1 = 4;    // explicitly add in terms for l1=4, l2>=4
                    int c  = i + l1 + 1;
                    for (int l2 = 4; l2 <= j - e - 2; ++l2) {
                        int s = l1 + l2;
                        int f = j - l2 - 1;
                        if (allow_paired(c, f) && Qx(i, d, e, s) > 0.0) {
                            float temp =
                            Qg(c, d, e, f) *
                            EXP(-(score_interior_asymmetry(l1, l2) + score_interior_mismatch(f, c, f + 1, c - 1)) / RT);
                            float p = temp / Qx(i, d, e, s) * Px(i, d, e, s);
                            Pg(c, d, e, f) += p;
                            Px(i, d, e, s) -= p;
                            assert(!std::isnan(p));
                            if (temp > Qx(i, d, e, s)) {
                                temp      = 0.0;
                                int l1min = 5;
                                int l2min = 4;
                                for (int c = i + l1min + 1; c <= d - 1; ++c) {
                                    int f = c - i + j - 3;
                                    if (j - f - 1 > l2min && f >= e + 1 && allow_paired(c, f)) {
                                        temp += Qg(c, d, e, f) * EXP(-score_interior(i, j, c, f, true) / RT);
                                    }
                                }
                                Qx(i, d, e, s) = temp;
                            } else {
                                Qx(i, d, e, s) -= temp;
                            }
                        }
                    }

                    if (d >= i + 7) {
                        int l2 = 4;    // explicitly add in terms of l1>=5, l2=4
                        int f  = j - l2 - 1;
                        for (int l1 = 5; l1 <= d - i - 2; ++l1) {
                            int s = l1 + l2;
                            int c = i + l1 + 1;
                            if (allow_paired(c, f) && Qx(i, d, e, s) > 0.0) {
                                float temp =
                                Qg(c, d, e, f) *
                                EXP(-(score_interior_asymmetry(l1, l2) + score_interior_mismatch(f, c, f + 1, c - 1)) / RT);
                                float p = temp / Qx(i, d, e, s) * Px(i, d, e, s);
                                Pg(c, d, e, f) += p;
                                Px(i, d, e, s) -= p;
                                assert(!std::isnan(p));
                                if (temp > Qx(i, d, e, s)) {
                                    temp      = 0.0;
                                    int l1min = 5;
                                    int l2min = 5;
                                    for (int c = i + l1min + 1; c <= d - 1; ++c) {
                                        int f = c - i + j - 4;
                                        if (j - f - 1 > l2min && f >= e + 1 && allow_paired(c, f)) {
                                            temp += Qg(c, d, e, f) * EXP(-score_interior(i, j, c, f, true) / RT);
                                        }
                                    }
                                    Qx(i, d, e, s) = temp;
                                } else {
                                    Qx(i, d, e, s) -= temp;
                                }
                            }
                        }
                    }
                }

                // Store paRTial values for Qx2 and Px2
                for (int s = 10; s <= l - 9; ++s) {
                    Qx2(i + 1, d, e, s - 2) = Qx(i, d, e, s) * EXP(-(score_loop(s + 2) - score_loop(s)) / RT);
                    Px2(i + 1, d, e, s - 2) = Px(i, d, e, s);
                }
            }
        }
    }
}

void Nupack::get_posterior(vector<float>& bp1, vector<float>& bp2, vector<int>& offset) const
{
    bp1.resize((N + 1) * (N + 2) / 2);
    bp2.resize((N + 1) * (N + 2) / 2);
    offset.resize(N + 1);
    for (int i = 0; i <= N; ++i) offset[i] = i * ((N + 1) + (N + 1) - i - 1) / 2;
    for (int i = 0; i != N - 1; ++i)
        for (int j = i + 1; j != N; ++j) {
            bp1[offset[i + 1] + (j + 1)] = Pb(i, j);
            bp2[offset[i + 1] + (j + 1)] = Pbg(i, j);
        }
}

void Nupack::get_posterior(vector<float>& bp, vector<int>& offset) const
{
    bp.resize((N + 1) * (N + 2) / 2);
    offset.resize(N + 1);
    for (int i = 0; i <= N; ++i) offset[i] = i * ((N + 1) + (N + 1) - i - 1) / 2;
    for (int i = 0; i != N - 1; ++i)
        for (int j = i + 1; j != N; ++j) bp[offset[i + 1] + (j + 1)] = Pb(i, j) + Pbg(i, j);
}

energy_t Nupack::score_hairpin(int i, int j) const
{
    energy_t e     = 0.0;
    bool     polyC = true;
    for (int k = i + 1; k < j; ++k) {
        if (seq[k] != BASE_C) {
            polyC = false;
            break;
        }
    }

    int size = j - i - 1;

    assert(size >= 3);
    assert(allow_paired(i, j));

    e += size <= 30 ? hairpin37[size - 1] : hairpin37[30 - 1] + loop_greater30 * LOG(size / 30.0);

    if (size == 3) {
        e += score_at_penalty(i, j);
        e += triloop37[seq[i] - 1][seq[i + 1] - 1][seq[i + 2] - 1][seq[j - 1] - 1][seq[j] - 1];
        if (polyC) e += polyC_penalty;
        if (seq[i + 1] == BASE_G && seq[i + 2] == BASE_G && seq[j - 1] == BASE_G) e += hairpin_GGG;
    } else if (size == 4) {
        e += tloop37[seq[i] - 1][seq[i + 1] - 1][seq[i + 2] - 1][seq[j - 2] - 1][seq[j - 1] - 1][seq[j] - 1];
        e += mismatch_hairpin37[seq[i + 1] - 1][seq[j - 1] - 1][pair_type(i, j)];
        if (polyC) e += polyC_slope * size + polyC_int;
    } else /*if (size>4)*/
    {
        e += mismatch_hairpin37[seq[i + 1] - 1][seq[j - 1] - 1][pair_type(i, j)];
        if (polyC) e += polyC_slope * size + polyC_int;
    }
    return e;
}

energy_t Nupack::score_loop(int l) const
{
    return l <= 30 ? interior37[l - 1] : interior37[30 - 1] + loop_greater30 * LOG(l / 30.0);
}

energy_t Nupack::score_interior(int i, int h, int m, int j, bool pk) const
{
    int      l1   = h - i - 1;
    int      l2   = j - m - 1;
    int      size = l1 + l2;
    energy_t e    = 0;

    // helix
    if (size == 0) {
        return stack37[pair_type(i, j)][pair_type(h, m)] * (pk ? pk_stack_span : 1.0);
    }

    // bulge
    else if (l1 == 0 || l2 == 0) {
        e += size <= 30 ? bulge37[size - 1] : bulge37[30 - 1] + loop_greater30 * LOG(size / 30.0);

        if (l1 + l2 == 1)    // single bulge...treat as a stacked region
        {
            e += stack37[pair_type(i, j)][pair_type(h, m)];
            e -= salt_correction;
        } else {
            e += score_at_penalty(i, j);
            e += score_at_penalty(h, m);
        }
    }

    // interior loop
    else if (l1 > 0 && l2 > 0) {
        int asymmetry = std::abs(l1 - l2);
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
            e += int11_37[pair_type(i, j)][pair_type(h, m)][seq[i + 1] - 1][seq[j - 1] - 1];
        else if (l1 == 2 && l2 == 2)
            e += int22_37[pair_type(i, j)][pair_type(h, m)][seq[i + 1] - 1][seq[j - 1] - 1][seq[i + 2] - 1][seq[j - 2] - 1];
        else if (l1 == 1 && l2 == 2)
            e += int21_37[pair_type(i, j)][seq[j - 2] - 1][seq[i + 1] - 1][pair_type(h, m)][seq[j - 1] - 1];
        else if (l1 == 2 && l2 == 1)
            e += int21_37[pair_type(m, h)][seq[i + 1] - 1][seq[j - 1] - 1][pair_type(j, i)][seq[i + 2] - 1];
        else {
            assert(!"error in tabulated interior loop");
            exit(1);
        }
    } else {
        assert(!"improperly classifed interior loop");
        exit(1);
    }
    return e * (pk ? pk_interior_span : 1.0);
}

energy_t Nupack::score_interior_mismatch(int i, int j, int k, int l) const
{
    return mismatch_interior37[seq[k] - 1][seq[l] - 1][pair_type(i, j)];
}

energy_t Nupack::score_interior_mismatch(int i, int j) const
{
    return mismatch_interior37[BASE_N][BASE_N][pair_type(i, j)];
}

energy_t Nupack::score_interior_asymmetry(int l1, int l2) const
{
    energy_t e         = 0.0;
    int      size      = l1 + l2;
    int      asymmetry = std::abs(l1 - l2);
    e += size <= 30 ? interior37[size - 1] : interior37[30 - 1] + loop_greater30 * LOG(size / 30.0);

    // asymmetry penalty
    e += std::min(max_asymmetry, asymmetry * asymmetry_penalty[std::min(4, std::min(l1, l2)) - 1]);

    return e;
}

energy_t Nupack::score_multiloop(bool pk) const { return pk ? multiloop_penalty_pk : multiloop_penalty; }

energy_t Nupack::score_multiloop_paired(int n, bool pk) const
{
    return (pk ? multiloop_paired_penalty_pk : multiloop_paired_penalty) * n;
}

energy_t Nupack::score_multiloop_unpaired(int n, bool pk) const
{
    return (pk ? multiloop_unpaired_penalty_pk : multiloop_unpaired_penalty) * n;
}

energy_t Nupack::score_at_penalty(int i, int j) const
{
    return pair_type(i, j) == PAIR_AU || pair_type(i, j) == PAIR_UA ? at_penalty : 0;
}

energy_t Nupack::score_dangle(int i, int j) const
{
    energy_t d5 = 0.0, d3 = 0.0;

    // if( DANGLETYPE != 2) {
    if ((j == i - 1) || (j == -1 && i > 0)) {
        return 0.0;
    }
    if ((j == -1 && i > 0) || (j == i - 1 && (i == 0 || j == N - 1))) return 0.0;
    if (j != N - 1) d3 = dangle3_37[3 - pair_type(j + 1)][seq[j] - 1];
    if (i != 0) d5 = dangle5_37[pair_type(i - 1)][seq[i] - 1];
    if (i == j && i != 0 && j != N - 1 /* && DANGLETYPE!=2 */)
        return std::min(d3, d5);
    else
        return d3 + d5;
}

inline energy_t Nupack::score_pk() const { return pk_penalty; }
inline energy_t Nupack::score_pk_multiloop() const { return pk_multiloop_penalty; }
inline energy_t Nupack::score_pk_pk() const { return pk_pk_penalty; }
inline energy_t Nupack::score_pk_paired(int n) const { return pk_paired_penalty * n; }
inline energy_t Nupack::score_pk_unpaired(int n) const { return pk_unpaired_penalty * n; }
inline energy_t Nupack::score_pk_band(int n) const { return pk_band_penalty * n; }
