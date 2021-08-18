#ifndef MOIP_H_
#define MOIP_H_

#define IL_STD

#include "SecondaryStructure.h"
#include "rna.h"
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <vector>

using std::vector;

typedef struct args_ {
						path           motif_file;
						std::mutex&    posInsertionSites_mutex;
						args_(path motif_file_, mutex& mutex_) : motif_file(motif_file_), posInsertionSites_mutex(mutex_) {}
					  } args_of_parallel_func;


class MOIP
{
	public:
	MOIP(void);
	MOIP(const RNA& rna, string source, string source_path, float theta, bool verbose);
	~MOIP(void);
	SecondaryStructure        	solve_objective(int o, double min, double max);
	SecondaryStructure        	solve_objective(int o);
	uint						get_n_candidates(void) const;
	uint                      	get_n_solutions(void) const;
	const SecondaryStructure& 	solution(uint i) const;
	void                      	search_between(double lambdaMin, double lambdaMax);
	bool                      	allowed_basepair(size_t u, size_t v) const;
	void                      	add_solution(const SecondaryStructure& s);
	void                      	remove_solution(uint i);
	void                      	forbid_solutions_between(double min, double max);
	IloEnv&                   	get_env(void);
	static char               	obj_function_nbr_;    // On what criteria do you want to insert motifs ?
	static uint               	obj_to_solve_;  // What objective do you prefer to solve in mono-objective portions of the algorithm ?
	static double             	precision_;   // decimals to keep in objective values, to avoid numerical issues. otherwise, solution with objective 5.0000000009 dominates solution with 5.0 =(
	static bool               	allow_pk_;      // Wether we forbid pseudoknots (false) or allow them (true)
	static uint               	max_sol_nbr_;  // Number of solutions to accept in the Pareto set before we give up the computation
	
	private:
	bool   						is_undominated_yet(const SecondaryStructure& s);
	void   						define_problem_constraints(string& source);
	size_t 						get_yuv_index(size_t u, size_t v) const;
	size_t 						get_Cpxi_index(size_t x_i, size_t i_on_j) const;
	IloNumExprArg& 				y(size_t u, size_t v);    // Direct reference to y^u_v in basepair_dv_
	IloNumExprArg& 				C(size_t x, size_t i);    // Direct reference to C_p^xi in insertion_dv_
	bool   						exists_vertical_outdated_labels(const SecondaryStructure& s) const;
	bool   						exists_horizontal_outdated_labels(const SecondaryStructure& s) const;
	void   						allowed_motifs_from_desc(args_of_parallel_func arg_struct);
	void   						allowed_motifs_from_rin(args_of_parallel_func arg_struct);
	void						allowed_motifs_from_json(args_of_parallel_func arg_struct, vector<pair<uint, char>> errors_id);
	
	bool verbose_;    // Should we print things ?

	// Elements of the problem
	RNA                        rna_;                // RNA object
	vector<Motif>              insertion_sites_;    // Potential Motif insertion sites
	vector<SecondaryStructure> pareto_;             // Vector of results

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

	vector<vector<uint>>   index_of_xij_;		         //Stores the indexes of the xij variables (BioKop)
};

inline uint                      MOIP::get_n_solutions(void) const { return pareto_.size(); }
inline uint                      MOIP::get_n_candidates(void) const { return insertion_sites_.size(); }
inline const SecondaryStructure& MOIP::solution(uint i) const { return pareto_[i]; }
inline IloNumExprArg&            MOIP::y(size_t u, size_t v) { return basepair_dv_[get_yuv_index(u, v)]; }
inline IloNumExprArg&            MOIP::C(size_t x, size_t i) { return insertion_dv_[get_Cpxi_index(x, i)]; }
inline SecondaryStructure        MOIP::solve_objective(int o) { return solve_objective(o, 0, rna_.get_RNA_length()); }
inline IloEnv&                   MOIP::get_env(void) { return env_; }

#endif    // MOIP_H_