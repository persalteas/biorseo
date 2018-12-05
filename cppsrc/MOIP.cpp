#include "MOIP.h"
#include "SecondaryStructure.h"
#include "rna.h"
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

uint MOIP::ncores = 0;

MOIP::MOIP(const RNA& rna, const vector<Motif>& insertionSites, float pthreshold)
: rna_(rna), insertion_sites_(insertionSites), theta_{pthreshold}
{
    basepair_dv_  = IloNumVarArray(env_);
    insertion_dv_ = IloNumVarArray(env_);

    // Add the y^u_v decision variables
    uint u, v, c = 0;
    index_of_yuv_ = vector<vector<size_t>>(rna_.get_RNA_length() - 6, vector<size_t>(0));
    for (u = 0; u < rna_.get_RNA_length() - 6; u++) {
        for (v = u + 4; v < rna_.get_RNA_length(); v++)    // A basepair is possible iff v > u+3
        {
            if (rna_.get_pij(u, v) > theta_) {
                index_of_yuv_[u].push_back(c);
                c++;
                char name[15];
                sprintf(name, "y%d,%d", u, v);
                basepair_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether u and v are paired
            } else {
                index_of_yuv_[u].push_back(rna_.get_RNA_length() + 1);
            }
        }
    }
    // Add the Cx,i,p decision variables
    index_of_first_components.reserve(insertionSites.size());
    index_of_Cxip_.reserve(insertionSites.size());
    size_t i = 0;
    for (const Motif m : insertionSites) {
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
            insertion_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether component i of
                                                                                // motif x is inserted at position p
        }
    }

    rna_.print_basepair_p_matrix(theta_);

    cout << c << " + " << i << " (yuv + Cpxi) decision variables are used !" << endl;
}

MOIP::~MOIP() { env_.end(); }



bool MOIP::is_undominated_yet(const SecondaryStructure& s)
{
    for (uint i = 0; i < pareto_.size(); i++) {
        if (pareto_[i] > s) return false;
    }
    return true;
}

SecondaryStructure MOIP::solve_objective(int o, double min, double max)
{
    IloModel model_ = IloModel(env_);
    cout << "Solving objective function " << o << "..." << endl;
    add_problem_constraints(model_);
    switch (o) {
    case 1: {
        // Add the motif objective function
        IloExpr obj1 = IloExpr(env_);
        for (uint i = 0; i < insertion_sites_.size(); i++) {
            IloNum n_compo_squared = IloNum(insertion_sites_[i].comp.size() * insertion_sites_[i].comp.size());
            obj1 += n_compo_squared * insertion_dv_[index_of_first_components[i]];
        }
        model_.add(IloMaximize(env_, obj1));
    } break;
    case 2: {
        // Add the likelihood objective function
        IloExpr obj2 = IloExpr(env_);
        for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
            for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
                if (allowed_basepair(u, v)) obj2 += (IloNum(rna_.get_pij(u, v)) * y(u, v));
            }
        }
        model_.add(IloMaximize(env_, obj2));
    }
    }
    IloCplex cplex_ = IloCplex(model_);
    if (!cplex_.solve()) {
        env_.error() << "\t>Failed to optimize LP." << endl;
        throw(-1);
    }
    IloNumArray basepair_values(env_);
    IloNumArray insertion_values(env_);
    env_.out() << endl << "Solution status = " << cplex_.getStatus() << endl;
    env_.out() << endl << "Objective value = " << cplex_.getObjValue() << endl;
    cplex_.getValues(basepair_values, basepair_dv_);
    cplex_.getValues(insertion_values, insertion_dv_);

    // Build a secondary Structure
    SecondaryStructure best_ss = SecondaryStructure(rna_);
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        // A constraint requires that all the components are inserted or none, so testing the first is enough:
        if (insertion_values[index_of_first_components[i]]) best_ss.insert_motif(insertion_sites_[i]);
    }
    cout << "\t>retrieveing motifs inserted in the result secondary structure..." << endl;
    for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
        for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
            if (allowed_basepair(u, v))
                if (basepair_values[get_yuv_index(u, v)]) best_ss.set_basepair(u, v);
        }
    }
    cout << "\t>retrieving basepairs of the result secondary structure..." << endl;
    return best_ss;
}

void MOIP::add_problem_constraints(const IloModel& model_)
{
    // ensure there only is 0 or 1 pairing by nucleotide:
    cout << "\t>ensuring there are 0 or 1 pairing by nucleotide..." << endl;
    uint u, v;
    uint n = rna_.get_RNA_length();
    for (u = 0; u < n - 6; u++) {
        IloExpr c1(env_);
        for (v = 0; v < u; v++)
            if (allowed_basepair(v, u)) c1 += y(v, u);
        for (v = u + 4; v < n; v++)
            if (allowed_basepair(u, v)) c1 += y(u, v);
        model_.add(c1 <= 1);
        // cout << "\t>It worked for base " << u << " : " << (c1 <= 1) << endl;
    }
    // forbid lonely basepairs
    cout << "\t>forbidding lonely basepairs..." << endl;
    for (u = 0; u < n - 6; u++) {
        IloExpr c2(env_);    // for the case where s[u] is paired to s[v], v>u
        c2 += 1;
        for (v = u; v < n; v++)
            if (allowed_basepair(u - 1, v)) c2 += y(u - 1, v);
        for (v = u + 1; v < n; v++)
            if (allowed_basepair(u, v)) c2 -= y(u, v);
        for (v = u + 2; v < n; v++)
            if (allowed_basepair(u + 1, v)) c2 += y(u + 1, v);
        model_.add(c2 >= 1);
        // cout << "\t>It worked for base " << u << " : " << (c2 >= 1) << endl;
    }
    for (v = 5; v < n; v++) {
        IloExpr c2p(env_);    // for the case where s[u] is paired to s[v], v<u
        c2p += 1;
        for (u = 1; u <= v - 2; u++)
            if (allowed_basepair(u, v - 1)) c2p += y(u, v - 1);
        for (u = 1; u <= v - 1; u++)
            if (allowed_basepair(u, v)) c2p -= y(u, v);
        for (u = 1; u <= v; u++)
            if (allowed_basepair(u, v + 1)) c2p += y(u, v + 1);
        model_.add(c2p >= 1);
        // cout << "\t>It worked for base " << u << " : " << (c2p >= 1) << endl;
    }
    // Forbid pairings inside every motif component if included
    cout << "\t>forbidding basepairs inside included motif's components..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif& x = insertion_sites_[i];
        for (size_t j = 0; j < x.comp.size(); j++) {
            Component& c = x.comp[j];
            IloExpr    c3(env_);
            IloNum     kxi = IloNum(c.k);
            c3 += kxi * C(i, j);
            for (u = c.pos.first; u < c.pos.second; u++) {
                for (v = 0; v < n; v++)
                    if (allowed_basepair(u, v)) c3 += y(u, v);
            }
            model_.add(c3 <= kxi);
        }
    }
    // Forbid component overlap
    cout << "\t>forbidding component overlap..." << endl;
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
        if (nterms) model_.add(c4 <= 1);
        // cout << "\t>It worked for base " << u << " : " << (c4 <= 1) << endl;
    }
    // Component completeness
    cout << "\t>ensuring that motives cannot be partially included..." << endl;
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
        // cout << "\t>It worked for motif " << i << " : " << (c5 == jm1 * C(i, 0)) << endl;
    }
    // Force basepairs between the end of a component and the beginning of the next
    cout << "\t>forcing basepairs between bounds of inserted components..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif&  x   = insertion_sites_[i];
        IloExpr c6p = IloExpr(env_);
        if (allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second))
            c6p += y(x.comp[0].pos.first, x.comp.back().pos.second);
        cout << "\t\t" << (C(i, 0) <= c6p) << " (" << allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second) << ")" << endl;
        model_.add(C(i, 0) <= c6p);
        if (x.comp.size() == 1)    // This constraint is for multi-component motives.
            continue;
        for (size_t j = 1; j < x.comp.size(); j++) {
            IloExpr c6 = IloExpr(env_);
            if (allowed_basepair(x.comp[j - 1].pos.second, x.comp[j].pos.first))
                c6 += y(x.comp[j - 1].pos.second, x.comp[j].pos.first);
            model_.add(C(i, j) <= c6);
            cout << "\t\t" << (C(i, j) <= c6) << " (" << allowed_basepair(x.comp[j - 1].pos.second, x.comp[j].pos.first)
                 << ")" << endl;
        }
    }
}

void MOIP::extend_pareto(double lambdaMin, double lambdaMax)
{
    if (lambdaMin >= lambdaMax) {
        cerr << "lambdaMax < lambdaMin, something went wrong !" << endl;
        exit(EXIT_FAILURE);
    }

    // for any SecondaryStructure in pareto_ such that the value of the second
    // objective is between lambdaMin and lambdaMax
    // a DIFF() constraint and a mirror constraint is added
    for (uint i = 0; i < pareto_.size(); i++) {
        // DIFF()
        if (
        (abs(pareto_[i].get_objective_score(2) - lambdaMin) < PRECISION or pareto_[i].get_objective_score(2) > lambdaMin) and
        (abs(pareto_[i].get_objective_score(2) - lambdaMax) < PRECISION or pareto_[i].get_objective_score(2) < lambdaMax)) {
            // ip.add_bj_ct(pareto_[i]);
        }
        // mirror
        // if (
        // (abs(pareto_[i].get_obj2M_() - lambdaMin) < PRECISION or
        // pareto_[i].get_obj2M_() > lambdaMin) and (abs(pareto_[i].get_obj2M_() -
        // lambdaMax) < PRECISION or pareto_[i].get_obj2M_() < lambdaMax)) {
        //     ip.add_bj_ct_M(pareto_[i]);
        // }
    }

    // SecondaryStructure s = solve_objective(1, lambdaMin, lambdaMax);

    // if (is_undominated_yet(s)) {
    //     // adding the SecondaryStructure s to the set pareto_
    //     add_solution(s);
    //     // run localPareto above the SecondaryStructure s
    //     extend_pareto(s.get_objective_score(2), lambdaMax);
    // }
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
    if (get_yuv_index(a, b) == rna_.get_RNA_length() + 1) return false;    // not allowed because proba < theta_
    return true;
}
