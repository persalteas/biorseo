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

using std::cerr, std::cout, std::endl;
using std::make_pair;
using std::vector;

uint MOIP::ncores = 0;

MOIP::MOIP(const RNA& rna, const vector<Motif>& insertionSites)
: rna_(rna), insertion_sites_(insertionSites), beta_(1.0), theta_{0.01}
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
    index_of_Cxip_ = vector<vector<size_t>>(0);
    index_of_Cxip_.reserve(insertionSites.size());
    size_t i = 0;
    for (const Motif m : insertionSites) {
        index_of_first_components.push_back(i);
        index_of_Cxip_.push_back(vector<size_t>(0));
        for (const Component cmp : m.comp) {
            index_of_Cxip_.back().push_back(i);
            if (cmp.k > 0) i++;
            char name[20];
            sprintf(
            name, "C%d,%d-%d", static_cast<int>(index_of_Cxip_.size() - 1),
            static_cast<int>(index_of_Cxip_.back().size() - 1), cmp.pos.first);
            insertion_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether component i of motif x is inserted at position p
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

void MOIP::solve_objective(int o, double min, double max)
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
        model_.add(IloMinimize(env_, obj1));
    } break;
    case 2: {
        // Add the likelihood objective function
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
    env_.out() << endl << "Basepair Values = " << basepair_values << endl;
    cplex_.getValues(insertion_values, basepair_dv_);
    env_.out() << endl << "Insertion Values = " << insertion_values << endl;

    // TODO : retrieve the secondary structure !!
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
    // To be continued ...
}

void MOIP::extend_pareto(double lambdaMin, double lambdaMax)
{
    if (lambdaMin >= lambdaMax) {
        cerr << "lambdaMax < lambdaMin, something's wrong !" << endl;
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