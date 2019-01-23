#include <iostream>
// extern "C" {
// #include <ViennaRNA/fold.h>
// #include <ViennaRNA/part_func.h>
// #include <ViennaRNA/utils/basic.h>
// }
// Import NUPACK energy computation
extern "C" {
#include "nupack/shared.h"
#include "nupack/thermo/core.h"
}
#include "rna.h"

extern DBL_TYPE* pairPrPbg;    // for pseudoknots
extern DBL_TYPE* pairPrPb;     // for pseudoknots
extern double    CUTOFF;

using namespace Eigen;
using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

RNA::RNA(string name, string seq)
: name_(name), seq_(seq), n_(seq.size()), pij_(MatrixXf::Zero(n_, n_))    // pair_map(Matrix<pair_t, 5, 5>::Constant(PAIR_OTHER)),
{
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
    }
    if (contains_T) cout << "\tWARNING: Thymines automatically replaced by uraciles.";
    if (unknown_chars.size() > 0) {
        cout << "\tWARNING: Unknown chars in input sequence ignored : ";
        for (char c : unknown_chars) cout << c << " ";
        cout << endl;
    }
    cout << "\t>sequence formatted" << endl;

    // // Compute using ViennaRNA
    // cout << "\t>computing pairing probabilities (ViennaRNA's algorithm)..." << endl;
    // const char      *cseq = seq.c_str();
    // vrna_fold_compound_t  *vc = vrna_fold_compound(cseq, NULL, VRNA_OPTION_PF);
    // char      *propensity = (char *)vrna_alloc(sizeof(char) * (strlen(cseq) + 1));
    // vrna_ep_t *ptr, *pair_probabilities = NULL;
    // float     en = vrna_pf_fold(cseq, propensity, &pair_probabilities);
    // printf("%s\n%s [ %6.2f ]\n", cseq, propensity, en);
    // for (ptr = pair_probabilities; ptr->i != 0; ptr++)
    // {
    //     cout << ptr->i << '\t' << ptr->j << '\t' << ptr->p << endl;
    //     if (ptr->j <= int(n_)) pij_(ptr->i,ptr->j) = ptr->p;
    // }
    // free(pair_probabilities);
    // free(propensity);
    // vrna_fold_compound_free(vc);

    // Compute using Nupack
    cout << "\t>computing pairing probabilities..." << endl;
    DBL_TYPE pf;
    int      length, tmpLength;
    int      i, j, q, r;
    char     seqChar[MAXSEQLENGTH];    // Complete sequence
    int      seqNum[MAXSEQLENGTH + 1];

    // Init parameters
    DANGLETYPE                         = 1;
    TEMP_K                             = 37.0 + ZERO_C_IN_KELVIN;
    SODIUM_CONC                        = 1.0;
    MAGNESIUM_CONC                     = 0.0;
    USE_LONG_HELIX_FOR_SALT_CORRECTION = 0;

    // Init seqChar
    strcpy(seqChar, seq.c_str());
    tmpLength = length = strlen(seqChar);
    convertSeq(seqChar, seqNum, tmpLength);
    int ns1, ns2;
    getSequenceLength(seqChar, &ns1);
    getSequenceLengthInt(seqNum, &ns2);

    pairPr    = (DBL_TYPE*)calloc((length + 1) * (length + 1), sizeof(DBL_TYPE));
    pairPrPbg = (DBL_TYPE*)calloc((length + 1) * (length + 1), sizeof(DBL_TYPE));
    pairPrPb  = (DBL_TYPE*)calloc((length + 1) * (length + 1), sizeof(DBL_TYPE));

    pf = pfuncFullWithSym(seqNum, 5, RNA37, DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, 1, 1, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
    printf("\t\t>Free energy: %.14Lf kcal/mol\n", -kB * TEMP_K * logl(pf));

    // cout << "\t\t>Base pair probabilities:" << endl;
    for (i = 0; i < length; i++) {
        for (j = i + 1; j < length; j++) {    // upper diagonal
            pij_(i, j) = pairPr[(length + 1) * i + j];
            // printf("\t\t%d %d\t%.4Le + %.4Le = %.4Le\t%s\n",i + 1,j + 1, pairPrPb[(length + 1) * i + j], pairPrPbg[(length
            // + 1) * i + j], pairPr[(length + 1) * i + j], pairPrPb[(length + 1) * i + j]+ pairPrPbg[(length + 1) * i + j] == pairPr[(length + 1) * i + j] ? "CHECK" : "ERROR");
        }
    }
    cout << endl << "\t\t>Fast checking..." << endl;
    vector<double> p_unpaired = vector<double>(n_, 0.0);
    for (i = 0; i < length; i++) {
        p_unpaired[i] = pairPr[(length + 1) * i + j];
        double sum    = 0.0;
        for (j = 0; j < length; j++) sum += pij_(i, j);
        printf("\t\t%d\tunpaired: %.4e\tpaired(pK+noPK): %.4e\tTotal: %f\n", i + 1, p_unpaired[i], sum, p_unpaired[i] + sum);
    }
    free(pairPr);
    free(pairPrPbg);
    free(pairPrPb);

    cout << "\t\t>pairing probabilities defined" << endl;
}

void RNA::print_basepair_p_matrix(float theta) const
{
    cout << endl << endl;
    cout << "\t=== -log10(p(i,j)) for each pair (i,j) of nucleotides: ===" << endl << endl;
    cout << "\t" << seq_ << endl;
    uint i = 0;
    for (uint u = 0; u < pij_.rows(); u++) {
        cout << "\t";
        for (uint v = 0; v < pij_.cols(); v++) {
            if (pij_(u, v) < 5e-10)
                cout << " ";
            else if (pij_(u, v) > theta)
                cout << "\033[0;32m" << int(-log10(pij_(u, v))) << "\033[0m";
            else
                cout << int(-log10(pij_(u, v)));
        }
        cout << seq_[i] << endl;
        i++;
    }
    cout << endl << "\t\033[0;32mgreen\033[0m basepairs are kept as decision variables." << endl << endl;
}

base_t RNA::base_type(char x) const
{
    if (x == 'a' or x == 'A') return BASE_A;
    if (x == 'c' or x == 'C') return BASE_C;
    if (x == 'g' or x == 'G') return BASE_G;
    if (x == 'u' or x == 'U') return BASE_U;
    return BASE_N;
}
