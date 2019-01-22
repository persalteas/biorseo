#include "SecondaryStructure.h"
#include <boost/format.hpp>

using std::abs;
using std::cout;
using std::endl;

static const double PRECISION(0.0001);

SecondaryStructure::SecondaryStructure(const RNA& rna)
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
    cout << '\t' << rna_.get_seq() << endl;
    cout << '\t' << to_string() << endl;
    for (const Motif& m : motif_info_) {
        uint i = 0;
        cout << '\t';
        for (auto c : m.comp) {
            while (i != c.pos.second + 1) {
                i++;
                if (i < c.pos.first + 1)
                    cout << " ";
                else
                    cout << '-';
            }
        }
        while (i < nBP_) {
            cout << " ";
            i++;
        }
        cout << "\t" << m.pos_string() << endl;
    }
    cout << endl;
}

void SecondaryStructure::sort(void)
{
    std::sort(basepairs_.begin(), basepairs_.end(), basepair_sorter);
    std::sort(motif_info_.begin(), motif_info_.end(), motif_sorter);
}

bool basepair_sorter(pair<uint, uint>& i, pair<uint, uint>& j)
{
    if (i.first < j.first) return true;
    if (i.first == j.first and i.second < j.second) return true;
    return false;
}

bool motif_sorter(Motif& m1, Motif& m2)
{
    if (m1.atlas_id.compare(m2.atlas_id) < 0) return true;
    return false;
}

bool operator>(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    double s11 = s1.get_objective_score(1);
    double s12 = s1.get_objective_score(2);
    double s21 = s2.get_objective_score(1);
    double s22 = s2.get_objective_score(2);

    bool obj1 = false, obj2 = false, strict1 = false, strict2 = false;

    if (s11 > s21) {
        strict1 = true;
        obj1    = true;
    } else if (s11 == s21) {
        obj1 = true;
    }
    if (s12 > s22) {
        strict2 = true;
        obj2    = true;
    } else if (s12 == s22) {
        obj2 = true;
    }

    if (obj1 && obj2 && (strict1 || strict2)) {
        return true;
    }

    return false;
}

bool operator<(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    double s11 = s1.get_objective_score(1);
    double s12 = s1.get_objective_score(2);
    double s21 = s2.get_objective_score(1);
    double s22 = s2.get_objective_score(2);

    bool obj1 = false, obj2 = false, strict1 = false, strict2 = false;

    if (s11 < s21) {
        strict1 = true;
        obj1    = true;
    } else if (s11 == s21) {
        obj1 = true;
    }
    if (s12 < s22) {
        strict2 = true;
        obj2    = true;
    } else if (s12 == s22) {
        obj2 = true;
    }

    if (obj1 && obj2 && (strict1 || strict2)) {
        return true;
    }
    return false;
}

bool operator>=(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    double s11 = s1.get_objective_score(1);
    double s12 = s1.get_objective_score(2);
    double s21 = s2.get_objective_score(1);
    double s22 = s2.get_objective_score(2);

    bool obj1 = false, obj2 = false, strict1 = false, strict2 = false;

    if (s11 > s21) {
        strict1 = true;
        obj1    = true;
    } else if (s11 == s21) {
        obj1 = true;
    }
    if (s12 > s22) {
        strict2 = true;
        obj2    = true;
    } else if (s12 == s22) {
        obj2 = true;
    }

    if (obj1 && obj2 && (strict1 || strict2) || (s11 == s21 && s12 == s22)) {
        return true;
    }

    return false;
}

bool operator<=(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    double s11 = s1.get_objective_score(1);
    double s12 = s1.get_objective_score(2);
    double s21 = s2.get_objective_score(1);
    double s22 = s2.get_objective_score(2);

    bool obj1 = false, obj2 = false, strict1 = false, strict2 = false;

    if (s11 < s21) {
        strict1 = true;
        obj1    = true;
    } else if (s11 == s21) {
        obj1 = true;
    }
    if (s12 < s22) {
        strict2 = true;
        obj2    = true;
    } else if (s12 == s22) {
        obj2 = true;
    }

    if (obj1 && obj2 && (strict1 || strict2) || (s11 == s21 && s12 == s22)) {
        return true;
    }
    return false;
}

bool operator==(const Component& c1, const Component& c2)
{
    if (c1.pos.first != c2.pos.first) return false;
    if (c1.pos.second != c2.pos.second) return false;
    if (c1.score != c2.score) return false;
    return true;
}

bool operator!=(const Component& c1, const Component& c2) { return not(c1 == c2); }

bool operator==(const Motif& m1, const Motif& m2)
{
    if (m1.atlas_id != m2.atlas_id) return false;
    if (m1.reversed != m2.reversed) return false;
    for (uint i = 0; i < m1.comp.size(); i++)
        if (m1.comp[i] != m2.comp[i]) return false;
    return true;
}

bool operator!=(const Motif& m1, const Motif& m2) { return not(m1 == m2); }

bool operator==(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    // Checks wether the secondary structures are exactly the same, including the inserted motifs.

    // fast checks to refute the equality
    if (s1.get_objective_score(1) != s2.get_objective_score(1)) return false;
    if (s1.get_objective_score(2) != s2.get_objective_score(2)) return false;
    if (s1.get_n_motifs() != s2.get_n_motifs()) return false;
    if (s1.get_n_bp() != s2.get_n_bp()) return false;

    // Deep checking
    for (uint i = 0; i < s1.get_n_bp(); i++)
        if (s1.basepairs_[i] != s2.basepairs_[i]) return false;
    for (uint i = 0; i < s1.get_n_motifs(); i++)
        if (s1.motif_info_[i] != s2.motif_info_[i]) return false;
    return true;
}