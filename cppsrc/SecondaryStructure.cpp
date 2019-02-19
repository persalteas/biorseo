#include "SecondaryStructure.h"
#include "MOIP.h"
#include <algorithm>
#include <boost/format.hpp>

using std::abs;
using std::cout;
using std::endl;

static const double PRECISION(0.0001);


SecondaryStructure::SecondaryStructure() {}


SecondaryStructure::SecondaryStructure(const RNA& rna)
: objective_scores_(vector<double>(2)), n_(rna.get_RNA_length()), nBP_(0), rna_(rna)
{
    is_empty_structure = false;
}

SecondaryStructure::SecondaryStructure(bool empty) : rna_(RNA()) { is_empty_structure = empty; }



string SecondaryStructure::to_DBN(void) const
{

    string                           res     = string(n_, '.');
    int                              pklevel = 0;
    bool                             crosses;
    char                             start, end;
    vector<vector<pair<uint, uint>>> noncrossingSG;    // Non crossing sets of edges, also called pseudoknot "levels"
    vector<uint>                     possible_sets;

    // get the non crossing subsets
    for (size_t i = 0; i < nBP_; i++) {
        possible_sets.clear();
        for (size_t j = 0; j < noncrossingSG.size(); j++) {
            // check if basepairs_[i] crosses the subset noncrossingSG[j]
            crosses = false;
            for (size_t k = 0; k < noncrossingSG[j].size(); k++) {
                if (
                ((basepairs_[i].first < noncrossingSG[j][k].first) and (basepairs_[i].second < noncrossingSG[j][k].second) and (basepairs_[i].second > noncrossingSG[j][k].first)) or
                ((basepairs_[i].first > noncrossingSG[j][k].first) and (basepairs_[i].second > noncrossingSG[j][k].second) and (noncrossingSG[j][k].second > basepairs_[i].first))) {
                    crosses = true;
                    break;
                }
            }
            // if not, consider possible to add basepairs_[i] to noncrossingSG[j]
            if (!crosses) possible_sets.push_back(j);
        }
        if (possible_sets.size()) {
            // add it to the largest possible subset (this is a heuristic, the real problem is a k-coloration NP hard problem)
            uint max = 0;
            uint pos = 0;
            for (size_t k = 0; k < possible_sets.size(); k++)
                if (max < noncrossingSG[possible_sets[k]].size()) {
                    max = noncrossingSG[possible_sets[k]].size();
                    pos = possible_sets[k];
                }
            noncrossingSG[pos].push_back(basepairs_[i]);
        } else {
            // If basepairs_[i] has not been inserted in any subset, create a new one
            noncrossingSG.push_back(vector<pair<uint, uint>>(1, basepairs_[i]));
        }
    }

    // get the sizes of the non crossing subsets
    vector<uint> SGsizes(noncrossingSG.size(), 0);
    for (size_t j = 0; j < noncrossingSG.size(); j++) SGsizes[j] = noncrossingSG[j].size();

    // Process the subsets from largest to thinest
    while (noncrossingSG.size()) {
        // Find the largest non-crossing subset
        uint j = std::distance(SGsizes.begin(), std::max_element(SGsizes.begin(), SGsizes.end()));
        // Apply basepairs in the output string
        for (uint i = 0; i < noncrossingSG[j].size(); i++) {
            switch (pklevel) {
            case 0: start = '(', end = ')'; break;
            case 1: start = '[', end = ']'; break;
            case 2: start = '{', end = '}'; break;
            case 3: start = '<', end = '>'; break;
            case 4: start = 'A', end = 'a'; break;
            case 5: start = 'B', end = 'b'; break;
            case 6: start = 'C', end = 'c'; break;
            case 7: start = 'D', end = 'd'; break;
            case 8: start = 'E', end = 'e'; break;
            case 9: start = 'F', end = 'f'; break;
            default: start = '|', end = '|';
            }
            res[noncrossingSG[j][i].first]  = start;
            res[noncrossingSG[j][i].second] = end;
        }
        pklevel++;
        // Remove the processed subset from the vector
        noncrossingSG.erase(std::remove(noncrossingSG.begin(), noncrossingSG.end(), noncrossingSG[j]), noncrossingSG.end());
        SGsizes.erase(SGsizes.begin() + j, SGsizes.begin() + j + 1);
    }
    return res;
}

string SecondaryStructure::to_string(void) const
{
    string s;
    s += to_DBN();
    for (const Motif& m : motif_info_) s += " + " + m.atlas_id;
    s += "\t" + boost::str(boost::format("%.7f") % objective_scores_[0]) + "\t" +
         boost::str(boost::format("%.7f") % objective_scores_[1]);
    return s;
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
        while (i < n_) {
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

    if (s11 - s21 > MOIP::precision_) {
        strict1 = true;
        obj1    = true;
    } else if (abs(s11 - s21) < MOIP::precision_) {
        obj1 = true;
    }
    if (s12 - s22 > MOIP::precision_) {
        strict2 = true;
        obj2    = true;
    } else if (abs(s12 - s22) < MOIP::precision_) {
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

    if (MOIP::precision_ < s21 - s11) {
        strict1 = true;
        obj1    = true;
    } else if (abs(s11 - s21) < MOIP::precision_) {
        obj1 = true;
    }
    if (MOIP::precision_ < s22 - s12) {
        strict2 = true;
        obj2    = true;
    } else if (abs(s12 - s22) < MOIP::precision_) {
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

    if (s11 - s21 > MOIP::precision_) {
        strict1 = true;
        obj1    = true;
    } else if (abs(s11 - s21) < MOIP::precision_) {
        obj1 = true;
    }
    if (s12 - s22 > MOIP::precision_) {
        strict2 = true;
        obj2    = true;
    } else if (abs(s12 - s22) < MOIP::precision_) {
        obj2 = true;
    }

    if ((obj1 && obj2 && (strict1 || strict2)) || ((abs(s11 - s21) < MOIP::precision_ && abs(s12 - s22) < MOIP::precision_))) {
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

    if (MOIP::precision_ < s21 - s11) {
        strict1 = true;
        obj1    = true;
    } else if (abs(s11 - s21) < MOIP::precision_) {
        obj1 = true;
    }
    if (MOIP::precision_ < s22 - s12) {
        strict2 = true;
        obj2    = true;
    } else if (abs(s12 - s22) < MOIP::precision_) {
        obj2 = true;
    }

    if ((obj1 && obj2 && (strict1 || strict2)) || ((abs(s11 - s21) < MOIP::precision_ && abs(s12 - s22) < MOIP::precision_))) {
        return true;
    }
    return false;
}

bool operator==(const Component& c1, const Component& c2)
{
    if (c1.pos.first != c2.pos.first) return false;
    if (c1.pos.second != c2.pos.second) return false;
    return true;
}

bool operator!=(const Component& c1, const Component& c2) { return not(c1 == c2); }



bool operator==(const Motif& m1, const Motif& m2)
{
    if (m1.atlas_id != m2.atlas_id) return false;
    if (m1.score != m2.score) return false;
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

bool operator!=(const SecondaryStructure& s1, const SecondaryStructure& s2)
{
    // Checks wether the secondary structures are different, including the inserted motifs.

    return not(s1 == s2);
}