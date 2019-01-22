#ifndef MOIP_H_
#define MOIP_H_

#define IL_STD

#include "SecondaryStructure.h"
#include "rna.h"
#include <ilcplex/ilocplex.h>

using std::vector;

const double PRECISION = 0.0001;


class MOIP
{
    public:
    static uint ncores;

    MOIP(const RNA& rna, const vector<Motif>& motifSites, float pthreshold);
    ~MOIP(void);
    SecondaryStructure        solve_objective(int o, double min, double max);
    SecondaryStructure        solve_objective(int o);
    uint                      get_n_solutions(void) const;
    const SecondaryStructure& solution(uint i) const;
    void                      extend_pareto(double lambdaMin, double lambdaMax);
    bool                      allowed_basepair(size_t u, size_t v) const;
    void                      add_solution(const SecondaryStructure& s);

    private:
    bool           is_undominated_yet(const SecondaryStructure& s);
    void           define_problem_constraints(void);
    size_t         get_yuv_index(size_t u, size_t v) const;
    size_t         get_Cpxi_index(size_t x_i, size_t i_on_j) const;
    IloNumExprArg& y(size_t u, size_t v);    // Direct reference to y^u_v in basepair_dv_
    IloNumExprArg& C(size_t x, size_t i);    // Direct reference to C_p^xi in insertion_dv_

    // Elements of the problem
    RNA                        rna_;                // RNA object
    vector<Motif>              insertion_sites_;    // Potential Motif insertion sites
    vector<SecondaryStructure> pareto_;             // Vector of results

    // Objectives related
    float theta_;    // theta parameter for the probability function

    // CPLEX objects
    IloEnv                 env_;                            // environment CPLEX object
    IloNumVarArray         basepair_dv_;                    // Decision variables
    IloNumVarArray         insertion_dv_;                   // Decision variables
    IloModel               model_;                      // Solver for objective 1
    IloExpr                obj1;                            // Objective function that counts inserted motifs
    IloExpr                obj2;                            // Objective function of expected accuracy
    vector<vector<size_t>> index_of_Cxip_;                  // Stores the indexes of the Cxip in insertion_dv_
    vector<size_t>         index_of_first_components;       // Stores the indexes of Cx1p in insertion_dv_
    vector<vector<size_t>> index_of_yuv_;                   // Stores the indexes of the y^u_v in basepair_dv_
};

inline uint                      MOIP::get_n_solutions(void) const { return pareto_.size(); }
inline const SecondaryStructure& MOIP::solution(uint i) const { return pareto_[i]; }
inline IloNumExprArg&            MOIP::y(size_t u, size_t v) { return basepair_dv_[get_yuv_index(u, v)]; }
inline IloNumExprArg&            MOIP::C(size_t x, size_t i) { return insertion_dv_[get_Cpxi_index(x, i)]; }
inline SecondaryStructure        MOIP::solve_objective(int o) { return solve_objective(o, 0, rna_.get_RNA_length()); }
#endif    // MOIP_H_