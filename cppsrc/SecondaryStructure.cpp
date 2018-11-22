#include "SecondaryStructure.h"
#include <boost/format.hpp>

static const double PRECISION(0.0001);

SecondaryStructure::SecondaryStructure(const vector<double>& scores, const vector<bool>& decision_variables, VII coord, int RNAlength)
: objective_scores_(scores), dv_(decision_variables), coord_(coord), n_(RNAlength)
{
}

string SecondaryStructure::to_DBN(void) const
{
    string res(n_, '.');
    for (size_t i = 0; i < n_; i++) {
        if (dv_[i]) {
            res[coord_[i].first]  = '(';
            res[coord_[i].second] = ')';
        }
    }
    return res;
}

string SecondaryStructure::to_string(void) const
{
    return to_DBN() + "\t" + boost::str(boost::format("%.6f") % objective_scores_[0]) + "\t" +
           boost::str(boost::format("%.6f") % objective_scores_[1]);
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
