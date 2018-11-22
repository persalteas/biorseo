#include "rna.h"
#include "nupack.h"
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::cout, std::cerr, std::endl;

RNA::RNA(string name, string seq)
{
    if (!check_seq(seq)) {
        cerr << "Unknown chars in input sequence. Please restrict to ACGTU." << endl;
        exit(EXIT_FAILURE);
    }
    name_ = name;
    seq_  = seq;
    format();
    n_ = seq_.size();
    cout << "\t>formatted sequence" << endl;

    /*define type_*/
    type_ = vector<vector<int>>(n_, vector<int>(n_));
    for (int i = 0; i < n_; i++) {
        for (int j = 0; j < n_; j++) {
            if (i < j) {
                std::stringstream ss;
                ss << seq_[i] << seq_[j];
                string str = ss.str();
                if (str.compare("AU") == 0) {
                    type_[i][j] = 1;
                } else if (str.compare("CG") == 0) {
                    type_[i][j] = 2;

                } else if (str.compare("GC") == 0) {
                    type_[i][j] = 3;
                } else if (str.compare("GU") == 0) {
                    type_[i][j] = 4;
                } else if (str.compare("UG") == 0) {
                    type_[i][j] = 5;
                } else if (str.compare("UA") == 0) {
                    type_[i][j] = 6;
                } else {
                    type_[i][j] = 0;
                }
            } else {
                type_[i][j] = 0;
            }
        }
    }
    nBP_ = type_.size();

    /*define coord_*/
    for (int i = 0; i < n_; i++) {
        for (int j = 0; j < n_; j++) {
            if (i < j and type_[i][j] > 0) {
                if (i != 0 and i != n_ and j != 0 and j != n_) {
                    if (type_[i - 1][j + 1] > 0 or type_[i + 1][j - 1] > 0) {
                        coord_.push_back(std::make_pair(i, j));
                    }
                } else if (i == 0 or j == n_) {
                    if (type_[i + 1][j - 1] > 0) {
                        coord_.push_back(std::make_pair(i, j));
                    }
                }
            }
        }
    }


    /*define pij_*/
    vector<float> bp;
    vector<int>   offset;
    Nupack        nu;
    // nu.load_parameters("rna1999.dG");
    nu.load_default_parameters();
    cout << "\t>default parameters loaded (Serra and Turner, 1995)" << endl;
    nu.load_sequence(seq_);
    cout << "\t>computing pairing probabilities..." << endl;
    try {
        nu.calculate_partition_function();
    } catch (std::exception& e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    nu.calculate_posterior();
    nu.get_posterior(bp, offset);

    pij_ = vector<vector<float>>(n_, vector<float>(n_));
    for (int i = 1; i <= n_; i++) {
        for (int j = 1; j <= n_; j++) {
            pij_[i - 1][j - 1] = bp[offset[i] + j];
        }
    }
    cout << "\t>pairing probabilities defined" << endl;
}

int RNA::find_coord(pair<int, int> p)
{
    vector<pair<int, int>>::iterator it = find(coord_.begin(), coord_.end(), p);
    int                              r  = -1;
    if (it != coord_.end()) r = distance(coord_.begin(), it);
    return r;
}

bool RNA::check_seq(string seq)    // Checks if the sequences only contains ACGUT.
{
    bool res = true;
    for (unsigned int i = 0; i < seq.size(); i++) {
        if (seq[i] != 'A' and seq[i] != 'U' and seq[i] != 'C' and seq[i] != 'G' and seq[i] != 'T') {
            res = false;
            break;
        }
    }
    return res;
}

void RNA::format()
{
    for (unsigned int i = 0; i < seq_.size(); i++) {
        seq_[i] = toupper(seq_[i]);
        if (seq_[i] == 'T') {
            seq_[i] = 'U';
            break;
        }
    }
}
