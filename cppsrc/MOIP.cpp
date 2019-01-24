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

MOIP::MOIP(const RNA& rna, const vector<Motif>& insertionSites, float pthreshold)
: rna_(rna), insertion_sites_(insertionSites), theta_{pthreshold}
{

    rna_.print_basepair_p_matrix(theta_);

    cout << "defining problem decision variables..." << endl;
    basepair_dv_  = IloNumVarArray(env_);
    insertion_dv_ = IloNumVarArray(env_);

    // Add the y^u_v decision variables
    cout << "\t>Legal basepairs : ";
    uint u, v, c = 0;
    index_of_yuv_ = vector<vector<size_t>>(rna_.get_RNA_length() - 6, vector<size_t>(0));
    for (u = 0; u < rna_.get_RNA_length() - 6; u++)
        for (v = u + 4; v < rna_.get_RNA_length(); v++)    // A basepair is possible iff v > u+3
            if (rna_.get_pij(u, v) > theta_) {
                cout << u << '-' << v << " ";
                index_of_yuv_[u].push_back(c);
                c++;
                char name[15];
                sprintf(name, "y%d,%d", u, v);
                basepair_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether u and v are paired
            } else {
                index_of_yuv_[u].push_back(rna_.get_RNA_length() + 1);
            }
    cout << endl;

    // Add the Cx,i,p decision variables
    cout << "\t>Candidate motif insertion sites : " << endl;
    index_of_first_components.reserve(insertionSites.size());
    index_of_Cxip_.reserve(insertionSites.size());
    size_t i = 0;
    for (const Motif m : insertionSites) {
        cout << "\t\t" << m.pos_string() << endl;
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

    cout << c << " + " << i << " (yuv + Cpxi) decision variables are used." << endl;

    // Adding the problem's constraints
    model_ = IloModel(env_);
    define_problem_constraints();

    // Define the motif objective function:
    obj1 = IloExpr(env_);
    for (uint i = 0; i < insertion_sites_.size(); i++) {
        IloNum n_compo_squared = IloNum(insertion_sites_[i].comp.size() * insertion_sites_[i].comp.size());
        obj1 += n_compo_squared * insertion_dv_[index_of_first_components[i]];
    }

    // Define the expected accuracy objective function:
    obj2 = IloExpr(env_);
    for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
        for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
            if (allowed_basepair(u, v)) obj2 += (IloNum(rna_.get_pij(u, v)) * y(u, v));
        }
    }

    // Add its objective function to each solver
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
    // Solves one of the objectives, under constraint that the other should be in [min, max]
    cout << "Solving objective function " << o << ", " << min << " <= Obj" << 3 - o << " <= " << max << "..." << endl;

    IloObjective obj;
    IloRange     bounds;

    switch (o) {
    case 1:
        obj    = IloMaximize(env_, obj1);
        bounds = IloRange(env_, min, obj2, max);
        break;
    case 2: obj = IloMaximize(env_, obj2); bounds = IloRange(env_, min, obj1, max);
    }
    model_.add(obj);
    model_.add(bounds);
    IloCplex cplex_ = IloCplex(model_);
    cplex_.setOut(env_.getNullStream());

    if (!cplex_.solve()) {
        env_.error() << "\t>Failed to optimize LP." << endl;
        throw(-1);
    }

    cout << "\t>Solution status: " << cplex_.getStatus() << ", with objective " << o << " value " << cplex_.getObjValue() << endl;
    IloNumArray basepair_values(env_);
    IloNumArray insertion_values(env_);
    cplex_.getValues(basepair_values, basepair_dv_);
    cplex_.getValues(insertion_values, insertion_dv_);

    cout << "\t>Building secondary structure..." << endl;

    // Build a secondary Structure
    SecondaryStructure best_ss = SecondaryStructure(rna_);
    // cout << "\t\t>retrieveing motifs inserted in the result secondary structure..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        // A constraint requires that all the components are inserted or none, so testing the first is enough:
        if (insertion_values[index_of_first_components[i]]) best_ss.insert_motif(insertion_sites_[i]);
    }
    // cout << "\t\t>retrieving basepairs of the result secondary structure..." << endl;
    for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
        for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
            if (allowed_basepair(u, v))
                if (basepair_values[get_yuv_index(u, v)]) best_ss.set_basepair(u, v);
        }
    }
    best_ss.sort();    // order the basepairs in the vector
    best_ss.set_objective_score(2, cplex_.getValue(obj2));
    best_ss.set_objective_score(1, cplex_.getValue(obj1));

    // Forbidding to find best_ss later : add a constraint to model_
    IloExpr c(env_);
    for (uint d = 0; d < insertion_dv_.getSize(); d++)
        if (cplex_.getValue(insertion_dv_[d]))
            c += IloNum(1) - insertion_dv_[d];
        else
            c += insertion_dv_[d];
    for (uint d = 0; d < basepair_dv_.getSize(); d++)
        if (cplex_.getValue(basepair_dv_[d]))
            c += IloNum(1) - basepair_dv_[d];
        else
            c += basepair_dv_[d];
    model_.add(c >= IloNum(1));

    // Removing the objective from the model_
    model_.remove(obj);
    model_.remove(bounds);

    return best_ss;
}

void MOIP::define_problem_constraints(void)
{
    // ensure there only is 0 or 1 pairing by nucleotide:
    cout << "\t>ensuring there are at most 1 pairing by nucleotide..." << endl;
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
            cout << "\t\t" << (c1 <= 1) << endl;
        }
    }
    // forbid lonely basepairs
    cout << "\t>forbidding lonely basepairs..." << endl;
    for (u = 0; u < n; u++) {
        IloExpr c2(env_);    // for the case where s[u] is paired to s[v], v>u
        count = 0;
        for (v = u; v < n; v++)
            if (allowed_basepair(u - 1, v)) c2 += y(u - 1, v);
        for (v = u + 1; v < n; v++)
            if (allowed_basepair(u, v)) {
                c2 -= y(u, v);
                count++;
            }
        for (v = u + 2; v < n; v++)
            if (allowed_basepair(u + 1, v)) c2 += y(u + 1, v);
        if (count) {
            model_.add(c2 >= 0);
            cout << "\t\t" << (c2 >= 0) << endl;
        }
    }
    for (v = 0; v < n; v++) {
        IloExpr c2p(env_);    // for the case where s[u] is paired to s[v], v<u
        count = 0;
        for (u = 0; u <= v - 2; u++)
            if (allowed_basepair(u, v - 1)) c2p += y(u, v - 1);
        for (u = 0; u <= v - 1; u++)
            if (allowed_basepair(u, v)) {
                c2p -= y(u, v);
                count++;
            }
        for (u = 0; u <= v; u++)
            if (allowed_basepair(u, v + 1)) c2p += y(u, v + 1);
        if (count) {
            model_.add(c2p >= 0);
            cout << "\t\t" << (c2p >= 0) << endl;
        }
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
            uint count = 0;
            for (u = c.pos.first + 1; u < c.pos.second - 1; u++) {
                for (v = 0; v < n; v++)
                    if (allowed_basepair(u, v)) {
                        c3 += y(u, v);
                        count++;
                    }
            }
            if (count) {
                model_.add(c3 <= kxi);
                cout << "\t\t";
                // cout << x.atlas_id << '-' << j << ": ";
                cout << (c3 <= kxi) << endl;
            }
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
        if (nterms) {
            model_.add(c4 <= 1);
            cout << "\t\t" << (c4 <= 1) << endl;
        }
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
        cout << "\t\t>motif " << i << " : " << (c5 == jm1 * C(i, 0)) << endl;
    }
    // Force basepairs between the end of a component and the beginning of the next
    cout << "\t>forcing basepairs between bounds of inserted components..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif&  x   = insertion_sites_[i];
        IloExpr c6p = IloExpr(env_);
        if (allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second))
            c6p += y(x.comp[0].pos.first, x.comp.back().pos.second);
        cout << "\t\t" << (C(i, 0) <= c6p) << "\t(" << x.comp[0].pos.first << ',' << x.comp.back().pos.second
             << (allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second) > 0 ? ") is allowed" : ") is not allowed")
             << endl;
        model_.add(C(i, 0) <= c6p);
        if (x.comp.size() == 1)    // This constraint is for multi-component motives.
            continue;
        for (size_t j = 0; j < x.comp.size()-1; j++) {
            IloExpr c6 = IloExpr(env_);
            if (allowed_basepair(x.comp[j].pos.second, x.comp[j+1].pos.first))
                c6 += y(x.comp[j].pos.second, x.comp[j+1].pos.first);
            model_.add(C(i, j) <= c6);
            cout << "\t\t" << (C(i, j) <= c6) << "\t(" << x.comp[j].pos.second << ',' << x.comp[j+1].pos.first
                 << (allowed_basepair(x.comp[j].pos.second, x.comp[j+1].pos.first) > 0 ? ") is allowed" : ") is not allowed")
                 << endl;
        }
    }
}

void MOIP::extend_pareto(double lambdaMin, double lambdaMax)
{
    SecondaryStructure s = solve_objective(1, lambdaMin, lambdaMax);
    cout << "\t>Done" << endl;
    if (is_undominated_yet(s)) {
        // adding the SecondaryStructure s to the set pareto_
        add_solution(s);
        // std::cin.ignore();
        // run localPareto above the SecondaryStructure s
        extend_pareto(s.get_objective_score(2), lambdaMax);
    }
}

void MOIP::add_solution(const SecondaryStructure& s)
{
    pareto_.push_back(s);
    for (uint i = 0; i < pareto_.size() - 1; i++)
        if (pareto_[i] < pareto_.back()) {
            // This should only happen in the case some structures have the same optimal Obj1 value.
            cout << "Removing structure from Pareto set : " << pareto_[i].to_string() << endl;
            pareto_.erase(pareto_.begin() + i);
        }
    cout << "Adding structure to Pareto set :     " << s.to_string() << endl;
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
