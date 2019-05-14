#include "MOIP.h"
#include <algorithm>
#include <boost/format.hpp>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using std::abs;
using std::cerr;
using std::cout;
using std::endl;
using std::make_pair;
using std::vector;

uint   MOIP::obj_function_nbr_ = 1;
uint   MOIP::obj_to_solve_     = 1;
double MOIP::precision_        = 1e-5;
bool   MOIP::allow_pk_         = true;

unsigned getNumConstraints(IloModel& m)
{
    unsigned           count = 0;
    IloModel::Iterator iter(m);
    while (iter.ok()) {
        if ((*iter).asConstraint().getImpl()) ++count;
        ++iter;
    }
    return count;
}

MOIP::MOIP() {}



MOIP::MOIP(const RNA& rna, const vector<Motif>& insertionSites, float theta, bool verbose)
: verbose_{verbose}, rna_(rna), insertion_sites_(insertionSites)
{

    if (verbose_) cout << "Summary of basepair probabilities:" << endl;
    if (verbose_) rna_.print_basepair_p_matrix(theta);

    if (verbose_) cout << "Defining problem decision variables..." << endl;
    basepair_dv_  = IloNumVarArray(env_);
    insertion_dv_ = IloNumVarArray(env_);

    // Add the y^u_v decision variables
    if (verbose_) cout << "\t>Legal basepairs : ";
    uint u, v, c = 0;
    index_of_yuv_ = vector<vector<size_t>>(rna_.get_RNA_length() - 6, vector<size_t>(0));
    for (u = 0; u < rna_.get_RNA_length() - 6; u++)
        for (v = u + 4; v < rna_.get_RNA_length(); v++)    // A basepair is possible iff v > u+3
            if (rna_.get_pij(u, v) > theta) {
                if (verbose_) cout << u << '-' << v << " ";
                index_of_yuv_[u].push_back(c);
                c++;
                char name[15];
                sprintf(name, "y%d,%d", u, v);
                basepair_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether u and v are paired
            } else {
                index_of_yuv_[u].push_back(rna_.get_RNA_length() * rna_.get_RNA_length() + 1);
            }
    if (verbose_) cout << endl;

    // Purge the list of insertion sites : remove those that are forbidden because bounding basepairs are forbidden
    if (verbose_) cout << "\t>Checking insertion sites..." << endl;
    vector<uint> to_remove;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif& x = insertion_sites_[i];
        if (verbose_) cout << "\t\tbasepair (" << x.comp[0].pos.first << ',' << x.comp.back().pos.second;
        if (allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second) > 0) {
            if (verbose_) cout << ") is allowed" << endl;
        } else if (!to_remove.size() or (to_remove.size() and to_remove.back() != i)) {
            if (verbose_) cout << ") is not allowed, removing motif " << i << " from candidates" << endl;
            to_remove.push_back(i);
            continue;
        } else {
            if (verbose_) cout << ") is not allowed (and motif has been previously forbidden)" << endl;
        }
        if (x.comp.size() == 1)    // This constraint is for multi-component motives.
            continue;
        for (size_t j = 0; j < x.comp.size() - 1; j++) {
            if (verbose_) cout << "\t\tbasepair (" << x.comp[j].pos.second << ',' << x.comp[j + 1].pos.first;
            if (allowed_basepair(x.comp[j].pos.second, x.comp[j + 1].pos.first) > 0) {
                if (verbose_) cout << ") is allowed" << endl;
            } else if (!to_remove.size() or (to_remove.size() and to_remove.back() != i)) {
                if (verbose_) cout << ") is not allowed, removing motif " << i << " from candidates" << endl;
                to_remove.push_back(i);
            } else {
                if (verbose_) cout << ") is not allowed (and motif has been previously forbidden)" << endl;
            }
        }
    }
    for (vector<uint>::reverse_iterator i = to_remove.rbegin(); i != to_remove.rend(); ++i)
        insertion_sites_.erase(insertion_sites_.begin() + (*i));

    // Add the Cx,i,p decision variables
    if (verbose_) cout << "\t>Allowed candidate insertion sites:" << endl;
    index_of_first_components.reserve(insertionSites.size());
    index_of_Cxip_.reserve(insertionSites.size());
    size_t i = 0;
    for (uint p = 0; p < insertion_sites_.size(); ++p) {
        const Motif& m = insertion_sites_[p];

        if (verbose_) cout << "\t\t>" << m.get_identifier() << '\t' << m.pos_string() << endl;
        index_of_first_components.push_back(i);
        index_of_Cxip_.push_back(vector<size_t>(0));
        for (const Component cmp : m.comp) {
            index_of_Cxip_.back().push_back(i);
            i++;
            char name[20];
            sprintf(
            name,
            "C%d,%d-%d",
            static_cast<int>(index_of_Cxip_.size() - 1),
            static_cast<int>(index_of_Cxip_.back().size() - 1),
            cmp.pos.first);
            insertion_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether component i of motif x is inserted at position p
        }
    }

    if (verbose_) cout << c << " + " << i << " (yuv + Cpxi) decision variables are used." << endl;

    // Adding the problem's constraints
    model_ = IloModel(env_);
    define_problem_constraints();
    if (verbose_) cout << "A total of " << getNumConstraints(model_) << " constraints are used." << endl;
    // if (getNumConstraints(model_) > 1500) {
    //     cerr << "\033[31m Quitting because too hard for me (too many constraints). Srry. \033[0m" << endl;
    //     exit(1);
    // }


    // if (getNumConstraints(model_) > 2000) {
    //     cerr << "\033[31mStopping 'cause too big for me...\033[0m" << endl;
    //     exit(-1);
    // }

    // Define the motif objective function:
    obj1 = IloExpr(env_);
    for (uint i = 0; i < insertion_sites_.size(); i++) {
        IloNum sum_k = 0;
        switch (obj_function_nbr_) {
        case 1:
            // RNA MoIP style
            for (const Component& c : insertion_sites_[i].comp) sum_k += c.k;
            obj1 += IloNum(sum_k * sum_k) * insertion_dv_[index_of_first_components[i]];
            break;

        case 2:
            // Weighted by the JAR3D score only:
            obj1 += IloNum(insertion_sites_[i].score_) * insertion_dv_[index_of_first_components[i]];
            break;

        case 3:
            // everything
            for (const Component& c : insertion_sites_[i].comp) sum_k += c.k;
            obj1 += IloNum(insertion_sites_[i].comp.size() * insertion_sites_[i].score_ / log2(sum_k)) *
                    insertion_dv_[index_of_first_components[i]];
            break;

        case 4:
            // everything but Jar3D score
            for (const Component& c : insertion_sites_[i].comp) sum_k += c.k;
            obj1 += IloNum(insertion_sites_[i].comp.size() / log2(sum_k)) * insertion_dv_[index_of_first_components[i]];
            break;
        }
    }

    // Define the expected accuracy objective function:
    obj2 = IloExpr(env_);
    for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
        for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
            if (allowed_basepair(u, v)) obj2 += (IloNum(rna_.get_pij(u, v)) * y(u, v));
        }
    }
}

MOIP::~MOIP() { env_.end(); }



bool MOIP::is_undominated_yet(const SecondaryStructure& s)
{
    for (SecondaryStructure& x : pareto_) {
        if (x > s) return false;
    }
    return true;
}

SecondaryStructure MOIP::solve_objective(int o, double min, double max)
{
    // Solves one of the objectives, under constraint that the other should be in [min, max]

    if (min > max) {
        // variable swap without a third, just because i want to look clever
        max = min + max;
        min = max - min;
        max = max - min;
    }

    // impose the bounds and the objective
    IloObjective obj;
    IloRange     bounds;
    switch (o) {
    case 1:
        obj    = IloMaximize(env_, obj1);
        bounds = IloRange(env_, min, obj2, max);
        break;
    case 2:
        obj    = IloMaximize(env_, obj2);
        bounds = IloRange(env_, min, obj1, max);
        break;
    }
    model_.add(obj);
    model_.add(bounds);

    IloCplex cplex_ = IloCplex(model_);
    cplex_.setOut(env_.getNullStream());
    // cplex_.exportModel("latestmodel.lp")

    if (!cplex_.solve()) {
        if (verbose_) cout << "\t>Failed to optimize LP: no more solutions to find." << endl;
        // Removing the objective from the model_
        model_.remove(obj);
        model_.remove(bounds);
        return SecondaryStructure(true);
    }

    if (verbose_)
        cout << "\t>Solution status: objective values (" << cplex_.getValue(obj1) << ", " << cplex_.getValue(obj2) << ')';

    // Build a secondary Structure
    SecondaryStructure best_ss = SecondaryStructure(rna_);
    // if (verbose_) cout << "\t\t>retrieveing motifs inserted in the result secondary structure..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++)
        // A constraint requires that all the components are inserted or none, so testing the first is enough:
        if (cplex_.getValue(insertion_dv_[index_of_first_components[i]]) > 0.5)
            best_ss.insert_motif(insertion_sites_[i]);

    // if (verbose_) cout << "\t\t>retrieving basepairs of the result secondary structure..." << endl;
    for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++)
        for (size_t v = u + 4; v < rna_.get_RNA_length(); v++)
            if (allowed_basepair(u, v))
                if (cplex_.getValue(y(u, v)) > 0.5) best_ss.set_basepair(u, v);

    best_ss.sort();    // order the basepairs in the vector
    best_ss.set_objective_score(2, cplex_.getValue(obj2));
    best_ss.set_objective_score(1, cplex_.getValue(obj1));

    // if (verbose_) cout << "\t\t>building the IP forbidding condition..." << endl;
    // Forbidding to find best_ss later
    IloExpr c(env_);
    for (uint d = 0; d < insertion_dv_.getSize(); d++)
        if (cplex_.getValue(insertion_dv_[d]) > 0.5)
            c += IloNum(1) - insertion_dv_[d];
        else
            c += insertion_dv_[d];
    for (uint d = 0; d < basepair_dv_.getSize(); d++)
        if (cplex_.getValue(basepair_dv_[d]) > 0.5)
            c += IloNum(1) - basepair_dv_[d];
        else
            c += basepair_dv_[d];
    model_.add(c >= IloNum(1));

    // exit
    model_.remove(bounds);
    model_.remove(obj);
    return best_ss;
}

void MOIP::define_problem_constraints(void)
{
    // Force basepairs between the end of a component and the beginning of the next
    if (verbose_) cout << "\t>forcing basepairs between bounds of inserted components..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif&  x   = insertion_sites_[i];
        IloExpr c6p = IloExpr(env_);
        if (allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second))
            c6p += y(x.comp[0].pos.first, x.comp.back().pos.second);
        if (verbose_) cout << "\t\t" << (C(i, 0) <= c6p) << endl;
        model_.add(C(i, 0) <= c6p);
        if (x.comp.size() == 1)    // This constraint is for multi-component motives.
            continue;
        for (size_t j = 0; j < x.comp.size() - 1; j++) {
            IloExpr c6 = IloExpr(env_);
            if (allowed_basepair(x.comp[j].pos.second, x.comp[j + 1].pos.first))
                c6 += y(x.comp[j].pos.second, x.comp[j + 1].pos.first);
            model_.add(C(i, j) <= c6);
            if (verbose_) cout << "\t\t" << (C(i, j) <= c6) << endl;
        }
    }

    // ensure there only is 0 or 1 pairing by nucleotide:
    if (verbose_) cout << "\t>ensuring there are at most 1 pairing by nucleotide..." << endl;
    uint u, v, count;
    uint n = rna_.get_RNA_length();
    for (u = 0; u < n; u++) {
        count = 0;
        IloExpr c1(env_);
        for (v = 0; v < u; v++)
            if (allowed_basepair(v, u)) {
                c1 += y(v, u);
                count++;
            }
        for (v = u + 4; v < n; v++)
            if (allowed_basepair(u, v)) {
                c1 += y(u, v);
                count++;
            }
        if (count > 1) {
            model_.add(c1 <= 1);
            if (verbose_) cout << "\t\t" << (c1 <= 1) << endl;
        }
    }
    // forbid lonely basepairs
    if (verbose_) cout << "\t>forbidding lonely basepairs..." << endl;
    for (u = 0; u < n - 5; u++)
        for (v = u + 4; v < n; v++) {
            if (allowed_basepair(u, v)) {
                IloExpr c2(env_);
                c2 += -y(u, v);
                if (allowed_basepair(u - 1, v + 1)) c2 += y(u - 1, v + 1);
                if (allowed_basepair(u + 1, v - 1)) c2 += y(u + 1, v - 1);
                model_.add(c2 >= 0);
                if (verbose_) cout << "\t\t" << (c2 >= 0) << endl;
            }
        }

    // Forbid pairings inside every motif component if included
    if (verbose_) cout << "\t>forbidding basepairs inside included motif's components..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif& x = insertion_sites_[i];
        for (size_t j = 0; j < x.comp.size(); j++) {
            Component& c = x.comp[j];
            IloExpr    c3(env_);
            IloNum     kxi = IloNum(c.k);
            c3 += (kxi - IloNum(2)) * C(i, j);
            uint count = 0;
            for (u = c.pos.first + 1; u < c.pos.second; u++) {
                for (v = 0; v < n; v++)
                    if (allowed_basepair(u, v)) {
                        c3 += y(u, v);
                        count++;
                    }
            }
            if (count > 0) {
                model_.add(c3 <= (kxi - IloNum(2)));
                if (verbose_) cout << "\t\t";
                if (verbose_) cout << x.get_identifier() << '-' << j << ": ";
                if (verbose_) cout << (c3 <= (kxi - IloNum(2))) << endl;
            }
        }
    }
    // Forbid component overlap
    if (verbose_) cout << "\t>forbidding component overlap..." << endl;
    for (u = 0; u < n; u++) {
        IloExpr c4(env_);
        uint    nterms = 0;
        for (size_t i = 0; i < insertion_sites_.size(); i++) {
            Motif& x = insertion_sites_[i];
            for (size_t j = 0; j < x.comp.size(); j++) {
                Component& c = x.comp[j];
                if (u >= c.pos.first and u <= c.pos.second) {    // Cxip contains u
                    c4 += C(i, j);
                    nterms++;
                }
            }
        }
        if (nterms > 1) {
            model_.add(c4 <= 1);
            if (verbose_) cout << "\t\t" << (c4 <= 1) << endl;
        }
    }
    // Component completeness
    if (verbose_) cout << "\t>ensuring that motives cannot be partially included..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif& x = insertion_sites_[i];
        if (x.comp.size() == 1)    // This constraint is for multi-component motives.
            continue;
        IloExpr c5(env_);
        IloNum  jm1 = IloNum(x.comp.size() - 1);
        for (size_t j = 1; j < x.comp.size(); j++) {
            c5 += C(i, j);
        }
        model_.add(c5 == jm1 * C(i, 0));
        if (verbose_) cout << "\t\t>motif " << i << " : " << (c5 == jm1 * C(i, 0)) << endl;
    }
    // Forbid pseudoknots
    if (!this->allow_pk_) {
        if (verbose_) cout << "\t>forbidding pseudoknots..." << endl;
        for (size_t u = 0; u < n - 6; u++)
            for (size_t v = u + 4; v < n - 1; v++)
                if (allowed_basepair(u, v))
                    for (size_t k = u + 1; k < v; ++k)
                        for (size_t l = v + 1; l < n; ++l)
                            if (allowed_basepair(k, l)) {
                                IloExpr c(env_);
                                c += y(u, v);
                                c += y(k, l);
                                model_.add(c <= 1);
                                if (verbose_) cout << "\t\t" << (c <= 1) << endl;
                            }
    }
}

void MOIP::search_between(double lambdaMin, double lambdaMax)
{
    SecondaryStructure s = solve_objective(obj_to_solve_, lambdaMin, lambdaMax);
    if (!s.is_empty_structure) {    // A solution has been found

        // if the solution is dominated, ignore it
        if (!is_undominated_yet(s)) {
            if (verbose_) cout << ", but structure is dominated." << endl;
            return;
        }

        // adding the SecondaryStructure s to the set pareto_
        if (verbose_) cout << ", not dominated." << endl;
        add_solution(s);

        // check if some labels should be updated on the vertical
        if (exists_vertical_outdated_labels(s))
            for (vector<SecondaryStructure>::iterator x = pareto_.end() - 2; x >= pareto_.begin(); x--)
                if (
                abs(x->get_objective_score(obj_to_solve_) - s.get_objective_score(obj_to_solve_)) < precision_ and
                precision_ < s.get_objective_score(3 - obj_to_solve_) - x->get_objective_score(3 - obj_to_solve_)) {
                    if (verbose_)
                        cout << "\t>removing structure from Pareto set, obj " << 3 - obj_to_solve_ << " = "
                             << x->get_objective_score(3 - obj_to_solve_) << endl;
                    pareto_.erase(x);
                }

        // search on top
        double min = s.get_objective_score(3 - obj_to_solve_) + precision_;
        double max = lambdaMax;
        if (verbose_)
            cout << std::setprecision(-log10(precision_) + 4) << "\nSolving objective function " << obj_to_solve_
                 << ", on top of " << s.get_objective_score(3 - obj_to_solve_) << ": Obj" << 3 - obj_to_solve_
                 << "  being in [" << std::setprecision(-log10(precision_) + 4) << min << ", "
                 << std::setprecision(-log10(precision_) + 4) << max << "]..." << endl;
        search_between(min, max);


        if (std::abs(max - min) - precision_ > precision_) {

            // search below
            min = lambdaMin;
            max = s.get_objective_score(3 - obj_to_solve_);
            if (verbose_)
                cout << std::setprecision(-log10(precision_) + 4) << "\nSolving objective function " << obj_to_solve_
                     << ", below (or eq. to) " << max << ": Obj" << 3 - obj_to_solve_ << "  being in ["
                     << std::setprecision(-log10(precision_) + 4) << min << ", "
                     << std::setprecision(-log10(precision_) + 4) << max << "]..." << endl;
            search_between(min, max);
        }

    } else {
        if (verbose_) cout << "\t>no solutions found." << endl;
    }
}

bool MOIP::exists_vertical_outdated_labels(const SecondaryStructure& s) const
{
    bool result = false;
    for (auto x : pareto_)
        if (x != s and abs(x.get_objective_score(obj_to_solve_) - s.get_objective_score(obj_to_solve_)) < precision_)
            result = true;
    if (result)
        for (auto x : pareto_)
            if (
            x != s and abs(x.get_objective_score(1) - s.get_objective_score(1)) < precision_ and
            abs(x.get_objective_score(2) - s.get_objective_score(2)) < precision_)
                result = false;
    return result;
}

bool MOIP::exists_horizontal_outdated_labels(const SecondaryStructure& s) const
{
    bool result = false;
    for (auto x : pareto_)
        if (x != s and abs(x.get_objective_score(3 - obj_to_solve_) - s.get_objective_score(3 - obj_to_solve_)) < precision_)
            result = true;
    if (result)
        for (auto x : pareto_)
            if (
            x != s and abs(x.get_objective_score(1) - s.get_objective_score(1)) < precision_ and
            abs(x.get_objective_score(2) - s.get_objective_score(2)) < precision_)
                result = false;
    return result;
}

void MOIP::add_solution(const SecondaryStructure& s)
{
    if (verbose_) cout << "\t>adding structure to Pareto set :\t" << s.to_string() << endl;
    pareto_.push_back(s);
    if (pareto_.size() > 500) {
        cerr << "\033[31m Quitting because combinatorial issues. \033[0m";
        exit(1);
    }
}

size_t MOIP::get_yuv_index(size_t u, size_t v) const
{
    size_t a, b;
    a = (u < v) ? u : v;
    b = (u > v) ? u : v;
    return index_of_yuv_[a][b - 4 - a];
}

size_t MOIP::get_Cpxi_index(size_t x_i, size_t i_on_j) const { return index_of_Cxip_[x_i][i_on_j]; }



bool MOIP::allowed_basepair(size_t u, size_t v) const
{
    size_t a, b;
    a = (v > u) ? u : v;
    b = (v > u) ? v : u;
    if (b - a < 4) return false;
    if (a >= rna_.get_RNA_length() - 6) return false;
    if (b >= rna_.get_RNA_length()) return false;
    if (get_yuv_index(a, b) == rna_.get_RNA_length() * rna_.get_RNA_length() + 1)
        return false;    // not allowed because proba < theta
    return true;
}

void MOIP::remove_solution(uint i) { pareto_.erase(pareto_.begin() + i); }