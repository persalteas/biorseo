#ifndef MOIP_H_
#define MOIP_H_

#define IL_STD

#include "SecondaryStructure.h"
#include "rna.h"
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>

using std::vector;

class MOIP
{
    public:
    MOIP(void);
    MOIP(const RNA& rna, const vector<Motif>& motifSites, uint nsets, float pthreshold, bool verbose);
    ~MOIP(void);
    SecondaryStructure        solve_objective(int o, double min, double max, const vector<IloConstraint>& F);
    SecondaryStructure        solve_objective(int o);
    uint                      get_n_solutions(void) const;
    const SecondaryStructure& solution(uint i) const;
    void                      search_between(double lambdaMin, double lambdaMax, const vector<IloConstraint>& F_);
    bool                      allowed_basepair(size_t u, size_t v) const;
    void                      add_solution(const SecondaryStructure& s);
    void                      remove_solution(uint i);
    vector<IloConstraint>     forbid_solutions_between(double min, double max);
    static uint obj_to_solve_;    // What objective do you prefer to solve in mono-objective portions of the algorithm ?
    static double precision_;    // decimals to keep in objective values, to avoid numerical issues. otherwise, solution with objective 5.0000000009 dominates solution with 5.0 =(
    static double epsilon_;

    private:
    bool   is_undominated_yet(const SecondaryStructure& s);
    void   define_problem_constraints(void);
    size_t get_yuv_index(size_t u, size_t v) const;
    size_t get_Cpxi_index(size_t x_i, size_t i_on_j) const;
    bool   exists_vertical_outdated_labels(const SecondaryStructure& s) const;
    bool   exists_horizontal_outdated_labels(const SecondaryStructure& s) const;

    IloNumExprArg& y(size_t u, size_t v);    // Direct reference to y^u_v in basepair_dv_
    IloNumExprArg& C(size_t x, size_t i);    // Direct reference to C_p^xi in insertion_dv_

    bool verbose_;    // Should we print things ?

    // Elements of the problem
    RNA                        rna_;                // RNA object
    vector<Motif>              insertion_sites_;    // Potential Motif insertion sites
    vector<SecondaryStructure> pareto_;             // Vector of results
    uint                       n_sets_;             // number of Pareto sets to return

    // Objectives related
    float theta_;    // theta parameter for the probability function

    // CPLEX objects
    IloEnv                 env_;                         // environment CPLEX object
    IloNumVarArray         basepair_dv_;                 // Decision variables
    IloNumVarArray         insertion_dv_;                // Decision variables
    IloModel               model_;                       // Solver for objective 1
    IloExpr                obj1;                         // Objective function that counts inserted motifs
    IloExpr                obj2;                         // Objective function of expected accuracy
    vector<vector<size_t>> index_of_Cxip_;               // Stores the indexes of the Cxip in insertion_dv_
    vector<size_t>         index_of_first_components;    // Stores the indexes of Cx1p in insertion_dv_
    vector<vector<size_t>> index_of_yuv_;                // Stores the indexes of the y^u_v in basepair_dv_
};

inline uint                      MOIP::get_n_solutions(void) const { return pareto_.size(); }
inline const SecondaryStructure& MOIP::solution(uint i) const { return pareto_[i]; }
inline IloNumExprArg&            MOIP::y(size_t u, size_t v) { return basepair_dv_[get_yuv_index(u, v)]; }
inline IloNumExprArg&            MOIP::C(size_t x, size_t i) { return insertion_dv_[get_Cpxi_index(x, i)]; }
inline SecondaryStructure        MOIP::solve_objective(int o)
{
    return solve_objective(o, 0, rna_.get_RNA_length(), vector<IloConstraint>());
}
#endif    // MOIP_H_