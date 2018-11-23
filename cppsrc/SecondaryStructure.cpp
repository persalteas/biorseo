#include "SecondaryStructure.h"
#include <boost/format.hpp>

using std::cout;
using std::endl;

static const double PRECISION(0.0001);

SecondaryStructure::SecondaryStructure(RNA& rna)
: objective_scores_(vector<double>(2)), n_(rna.get_RNA_length()), nBP_(0), rna_(rna)
{
}

string SecondaryStructure::to_DBN(void) const
{
    string res = string(n_, '.');
    for (size_t i = 0; i < nBP_; i++) {
        res[basepairs_[i].first]  = '(';
        res[basepairs_[i].second] = ')';
    }
    return res;
}

string SecondaryStructure::to_string(void) const
{
    return to_DBN() + "\t" + boost::str(boost::format("%.6f") % objective_scores_[0]) + "\t" +
           boost::str(boost::format("%.6f") % objective_scores_[1]);
}

void SecondaryStructure::set_basepair(uint i, uint j)
{
    nBP_++;
    pair<uint, uint> bp;
    bp.first  = i;
    bp.second = j;
    basepairs_.push_back(bp);
}

void SecondaryStructure::insert_motif(const Motif& m) { motif_info_.push_back(m); }

void SecondaryStructure::print(void) const
{
    cout << endl;
    cout << '\t' << to_DBN() << endl;
    for (const Motif& m : motif_info_) {
        uint i = 0;
        cout << '\t';
        for (auto c : m.comp) {
            while (i != c.pos.second) {
                if (i < c.pos.first)
                    cout << " ";
                else
                    cout << '-';
                i++;
            }
        }
        while (i < rna_.get_RNA_length()) {
            cout << " ";
            i++;
        }
        cout << "\t" << m.atlas_id << endl;
    }
    cout << endl;
}

bool operator>(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    double s11 = s1.get_objective_score(0);
    double s12 = s1.get_objective_score(1);
    double s21 = s2.get_objective_score(0);
    double s22 = s2.get_objective_score(1);

    bool obj1 = false, obj2 = false, strict1 = false, strict2 = false;

    if (s11 > s21) {
        strict1 = true;
        obj1    = true;
    } else if (abs(s11 - s21) < PRECISION) {
        obj1 = true;
    }
    if (s12 > s22) {
        strict2 = true;
        obj2    = true;
    } else if (abs(s12 - s22) < PRECISION) {
        obj2 = true;
    }

    if (obj1 && obj2 && (strict1 || strict2)) {
        return true;
    }

    return false;
}

bool operator<(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    double s11 = s1.get_objective_score(0);
    double s12 = s1.get_objective_score(1);
    double s21 = s2.get_objective_score(0);
    double s22 = s2.get_objective_score(1);

    bool obj1 = false, obj2 = false, strict1 = false, strict2 = false;

    if (s11 < s21) {
        strict1 = true;
        obj1    = true;
    } else if (abs(s11 - s21) < PRECISION) {
        obj1 = true;
    }
    if (s12 < s22) {
        strict2 = true;
        obj2    = true;
    } else if (abs(s12 - s22) < PRECISION) {
        obj2 = true;
    }

    if (obj1 && obj2 && (strict1 || strict2)) {
        return true;
    }
    return false;
}
