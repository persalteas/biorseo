#include "MOIP.h"
#include "Pool.h"
#include "Motif.h"
#include <algorithm>
#include <string> 
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <json.hpp>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

using namespace boost::filesystem;
using namespace std;
using json = nlohmann::json;

char   MOIP::obj_function_nbr_ = 'A';
char   MOIP::obj_function2_nbr_ = 'b';
uint   MOIP::obj_to_solve_     = 1;
double MOIP::precision_        = 1e-5;
bool   MOIP::allow_pk_         = true;
uint   MOIP::max_sol_nbr_      = 500;


struct recursive_directory_range {
    typedef recursive_directory_iterator iterator;
    recursive_directory_range(path p) : p_(p) {}

    iterator begin() { return recursive_directory_iterator(p_); }
    iterator end() { return recursive_directory_iterator(); }

    path p_;
};


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

MOIP::MOIP(const RNA& rna, string source, string source_path, float theta, bool verbose) : verbose_{verbose}, rna_(rna) 
{
    if (obj_function2_nbr_ != 'c' and !exists(source_path))
    {
        cerr << "ERR: Hmh, i can't find this: " << source_path << endl;
        exit(EXIT_FAILURE);
    }

    if (verbose_) cout << "Summary of basepair probabilities:" << endl;
    if (verbose_) rna_.print_basepair_p_matrix(theta);

    if (verbose_) cout << "Defining problem decision variables..." << endl;
    basepair_dv_  = IloNumVarArray(env_);
    insertion_dv_ = IloNumVarArray(env_);
    stacks_dv_ = IloNumVarArray(env_);

    // Add the y^u_v decision variables
    if (verbose_) cout << "\t> Legal basepairs : ";
    uint u, v, cy = 0;
    index_of_yuv_ = vector<vector<size_t>>(rna_.get_RNA_length() - 6, vector<size_t>(0));
    for (u = 0; u < rna_.get_RNA_length() - 6; u++)
        for (v = u + 4; v < rna_.get_RNA_length(); v++)    // A basepair is possible if v > u+3
            if (rna_.get_pij(u, v) > theta) {
                if (verbose_) cout << u << '-' << v << " ";
                index_of_yuv_[u].push_back(cy);
                cy++;
                char name[15];
                sprintf(name, "y%d,%d", u, v);
                basepair_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether u and v are paired
            } else {
                index_of_yuv_[u].push_back(rna_.get_RNA_length() * rna_.get_RNA_length() + 1);
            }
    if (verbose_) cout << endl;

    // Add the x_i,j decision variables
    if (verbose_) cout << "\t> The possible stacks of two base pairs (i,j) and (i+1,j-1) : ";
    uint cx = 0;
    index_of_xij_ = vector<vector<size_t>>(rna_.get_RNA_length() - 6, vector<size_t>(0));
    for (u = 0; u < rna_.get_RNA_length() - 6; u++)
        for (v = u + 4; v < rna_.get_RNA_length(); v++)    // A basepair is possible if v > u+3
            if (rna_.get_pij(u, v) > theta and rna_.get_pij(u + 1, v - 1) > theta) { // or u-1 v+1 ??
                if (verbose_) cout << u << '-' << v << " ";
                index_of_xij_[u].push_back(cx);
                cx++;
                char name[15];
                sprintf(name, "x%d,%d", u, v);
                stacks_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether (u,v) and (u+1,v-1) are a stack
            } else {
                index_of_xij_[u].push_back(rna_.get_RNA_length() * rna_.get_RNA_length() + 1);
            }
    if (verbose_) cout << endl;

    // Look for insertions sites, then create the appropriate Cxip variables
    insertion_sites_ = vector<Motif>();
    if (verbose_) cout << "\t> Looking for insertion sites..." << endl;

    if (source == "csvfile")
    {
        std::ifstream motifs;
        string        line;

        motifs = std::ifstream(source_path);
        getline(motifs, line);    // skip header
        while (getline(motifs, line))
        {

            Motif this_motif = Motif(line);
            bool to_keep = true;

            if (!(allowed_basepair(this_motif.comp[0].pos.first, this_motif.comp.back().pos.second)))
                // first nucleotide of first component and last nucleotide of last component cannot be paired,
                // so ignore this motif.
                continue;
            else if (this_motif.comp.size() != 1)
            {
                // Check that for every component, the last position can be paired to the next component's first
                for (size_t j = 0; j < this_motif.comp.size() - 1; j++)
                    if ( !(allowed_basepair(this_motif.comp[j].pos.second, this_motif.comp[j+1].pos.first)))
                    {
                        to_keep = false;
                        j = this_motif.comp.size(); // to exit the for loop()
                    }
                if (!to_keep) continue;
            }
            
            insertion_sites_.push_back(this_motif);
        }
    }
    else if (source == "descfolder") 
    {
        mutex         posInsertionSites_access;
        Pool          pool;
        int           errors   = 0;
        int           accepted = 0;
        int           inserted = 0;
        int           num_threads = thread::hardware_concurrency() - 1;
        vector<thread> thread_pool;

        for (int i = 0; i < num_threads; i++) 
            thread_pool.push_back(thread(&Pool::infinite_loop_func, &pool));

        // Read every .desc file and add it to the queue (iff valid)
        char error;
        for (auto it : recursive_directory_range(source_path))
        {    
            
            if ((error = Motif::is_valid_DESC(it.path().string()))) // Returns error if DESC file is incorrect
            {
                if (verbose)
                {
                    cerr << "\t> Ignoring motif " << it.path().stem();
                    switch (error)
                    {
                        case '-': cerr << ", some nucleotides have a negative number..."; break;
                        case 'l': cerr << ", hairpin (terminal) loops must be at least of size 3 !"; break;
                        case 'b': cerr << ", backbone link between non-consecutive residues ?"; break;
                        default:  cerr << ", use of an unknown nucleotide " << error;
                    }
                    cerr << endl;
                }
                errors++;
                continue;
            }
            accepted++;
            if (is_desc_insertible(it.path().string(), rna_.get_seq()))
            {
                file_and_mutex args(it.path(), posInsertionSites_access);
                inserted++;
                pool.push(bind(&MOIP::allowed_motifs_from_desc, this, args)); // & is necessary to get the pointer to a member function
            }
        }
        pool.done();

        for (unsigned int i = 0; i < thread_pool.size(); i++)
            thread_pool.at(i).join();

        if (verbose){
            cout << "\t> " << inserted << " candidate motifs on " << accepted + errors << " (" << errors << " ignored motifs), " << endl;
            cout << "\t  " << insertion_sites_.size() << " insertion sites kept after applying probability threshold of " << theta << endl;
        }
    }
    else if (source == "rinfolder")
    {
        mutex         posInsertionSites_access;
        Pool          pool;
        size_t        accepted = 0;
        size_t        errors   = 0;
        int           num_threads = thread::hardware_concurrency() - 1;
        vector<thread> thread_pool;

        for (int i = 0; i < num_threads; i++) 
            thread_pool.push_back(thread(&Pool::infinite_loop_func, &pool));

        // Read every RIN file and add it to the queue (iff valid)
        char error;
        for (auto it : recursive_directory_range(source_path))
        {
            if ((error = Motif::is_valid_RIN(it.path().string()))) // Returns error if RIN file is incorrect
            {
                if (verbose)
                {
                    cerr << "\t> Ignoring RIN " << it.path().stem();
                    switch (error)
                    {
                        case 'l': cerr << ", too short to be considered."; break;
                        case 'x': cerr << ", because not constraining the secondary structure."; break;
                        default: cerr << ", unknown reason";
                    }
                    cerr << endl;
                }
                errors++;
                continue;
            }
            accepted++;
            file_and_mutex args(it.path(), posInsertionSites_access);
            pool.push(bind(&MOIP::allowed_motifs_from_rin, this, args)); // & is necessary to get the pointer to a member function
        }
        pool.done();

        for (unsigned int i = 0; i < thread_pool.size(); i++)
            thread_pool.at(i).join();

        if (verbose){
            cout << "\t> " <<  insertion_sites_.size() << " candidate RINs on " << accepted + errors << " (" << errors << " ignored motifs), after applying probability threshold of " << theta << endl;
        }
    }
    else if (source == "jsonfolder")
    {
        mutex   posInsertionSites_access;
        Pool    pool;
        size_t  accepted = 0;
        size_t  errors   = 0;
        int     num_threads = thread::hardware_concurrency() - 1;
        vector<thread> thread_pool;
        std::ifstream  motif;

        // prepare a pool of threads to process motifs in parallel
        for (int i = 0; i < num_threads; i++) 
            thread_pool.push_back(thread(&Pool::infinite_loop_func, &pool));

        // Read every JSON entry and add it to the queue (iff valid)
        motif = std::ifstream(source_path);
        json js = json::parse(motif);
        char error;
        for(auto it = js.begin(); it != js.end(); ++it) {
            if ((error = Motif::is_valid_JSON(it))) // Returns error if JSON entry is incorrect
            {
                if (verbose)
                {
                    cerr << "\t> Ignoring motif " << it.key(); // field key, i.e. motif identifier
                    switch (error)
                    {
                        case 'l': cerr << ", too short to be considered."; break;
                        case 'x': cerr << ", sequence and secondary structure are of different size."; break;
                        case 'e' : cerr << ", sequence is empty."; break;
                        case 'f' : cerr << ", 2D is empty."; break;
                        case 'r' : cerr << ", RNA sequence not recognized, only use ACGTU."; break;
                        case 'n' : cerr << ", brackets are not balanced."; break;
                        default: cerr << ", unknown reason";
                    }
                    cerr << endl;
                }
                errors++;
                continue;
            }
            accepted++;

            // add a new job to the pool, to run allowed_motifs_from_json(args)
            motif_and_mutex args(it, posInsertionSites_access);
            pool.push(bind(&MOIP::allowed_motifs_from_json, this, args)); // & is necessary to get the pointer to a member function
        }
        pool.done(); // we won't add new jobs

        // Wait for jobs to complete
        for (unsigned int i = 0; i < thread_pool.size(); i++)
            thread_pool.at(i).join();

        if (verbose_) cout << "\t> " <<  insertion_sites_.size() << " candidate motifs on " << accepted + errors << " (" << errors << " ignored motifs), after applying probability threshold of " << theta << endl;
    }
    else if (obj_function2_nbr_ != 'c')
    {
        cout << "Err: Unknown module source." << endl;
    }

    // Add the Cx,i,p decision variables
    if (verbose_) cout << "\t> Allowed candidate insertion sites:" << endl;
    index_of_first_components.reserve(insertion_sites_.size()); // to remember the place of first components in insertion_dv_
    index_of_Cxip_.reserve(insertion_sites_.size());  // One vector per insertion_site/module, these vectors containing indexes of their components's dv in insertion_dv_.
    size_t i = 0;
    for (uint p = 0; p < insertion_sites_.size(); ++p) {
        const Motif& m = insertion_sites_[p];

        if (verbose_) cout << "\t\t> " << m.get_identifier() << '\t' << m.pos_string() << '\t' << m.sec_struct() << endl;
        index_of_first_components.push_back(i);
        index_of_Cxip_.push_back(vector<size_t>(0)); // A vector of size 0 (empty)

        for (const Component& cmp : m.comp) {
            index_of_Cxip_.back().push_back(i); // Add i to the current module vector
            i++;
            char name[20];
            sprintf(
                name,
                "C%d,%d[%d,%d]",
                static_cast<int>(index_of_Cxip_.size() - 1), // The number of motifs to date
                static_cast<int>(index_of_Cxip_.back().size() - 1), // The number of components in the last motif to date
                cmp.pos.first,
                cmp.pos.second
            );
            insertion_dv_.add(IloNumVar(env_, 0, 1, IloNumVar::Bool, name));    // A boolean whether component i of motif x is inserted at position p, named 'name'
        }
    }

    if (verbose_) cout << "\t> " << cy << " + " << cx << " + " << i << " (yuv + xuv + Cpxi) decision variables are used." << endl;

    // Adding the problem's constraints
    model_ = IloModel(env_);
    define_problem_constraints(source);
    if (verbose_) cout << "A total of " << getNumConstraints(model_) << " constraints are used." << endl;
    
    // if (getNumConstraints(model_) > 1500) {
    //     cerr << "\033[31m Quitting because too hard for me (too many constraints). Srry. \033[0m" << endl;
    //     exit(1);
    // }


    // Define the objective functions
    obj1 = IloExpr(env_);
    if (obj_function2_nbr_ != 'c') {
        // Define the motif objective function:
        for (uint i = 0; i < insertion_sites_.size(); i++) {
            IloNum sum_k = 0;

            switch (obj_function_nbr_) {
            case 'A':
                // RNA MoIP style
                for (const Component& c : insertion_sites_[i].comp) sum_k += c.k;
                obj1 += IloNum(sum_k * sum_k) * insertion_dv_[index_of_first_components[i]];
                break;
                
            case 'B':
                // everything but the Jar3D/Bayespairing score
                for (const Component& c : insertion_sites_[i].comp) sum_k += c.k;
                obj1 += IloNum(insertion_sites_[i].comp.size() / log2(sum_k)) * insertion_dv_[index_of_first_components[i]];
                break;   

            case 'C':
                // Weighted by the JAR3D or BayesPairing score only:
                obj1 += IloNum(insertion_sites_[i].score_) * insertion_dv_[index_of_first_components[i]];
                break;

            case 'D':
                // everything
                for (const Component& c : insertion_sites_[i].comp) sum_k += c.k;
                obj1 += IloNum(insertion_sites_[i].comp.size() * insertion_sites_[i].score_ / log2(sum_k)) *
                        insertion_dv_[index_of_first_components[i]];
                break;
            }
        }
    } else {
        // user passed both --mea and --mfe, this is Biokop mode, obj1 will be MEA
        if (verbose_) cout << "Running in Biokop mode: MEA versus MFE, ignoring motifs." << endl;
        for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
            for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
                if (allowed_basepair(u, v)) obj1 += (IloNum(rna_.get_pij(u, v)) * y(u, v));
            }
        }
        // set obj_function2_nbr_ to a so that obj2 is set to MFE just here after.
        obj_function2_nbr_ = 'a';
    }
    
    //Stacking energy parameter matrix
    double energy[7][7] = {
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 1.1, 2.1, 2.2, 1.4, 0.9, 0.6},
                {0.0, 2.1, 2.4, 3.3, 2.1, 2.1, 1.4},
                {0.0, 2.2, 3.3, 3.4, 2.5, 2.4, 1.5},
                {0.0, 1.4, 2.1, 2.5, 1.3, 1.3, 0.5},
                {0.0, 0.9, 2.1, 2.4, 1.3, 1.3, 1.0},
                {0.0, 0.6, 1.4, 1.5, 0.5, 1.0, 0.3}
            }; 

    obj2 = IloExpr(env_);
    switch (obj_function2_nbr_) { 
        case 'a':
            // Define the MFE (Minimum Free Energy):
            for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
                for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
                    if (get_xij_index(u, v) != rna_.get_RNA_length() * rna_.get_RNA_length() + 1) {
                        uint type1 = rna_.get_type()[u][v];
                        uint type2 = rna_.get_type()[u + 1][v - 1];
                        obj2 += IloNum(energy[type1][type2]) * x(u, v);
                    }
                }
            }
            break;
        case 'b':
            // Define the expected accuracy objective function:
            //MEA:
            for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++) {
                for (size_t v = u + 4; v < rna_.get_RNA_length(); v++) {
                    if (allowed_basepair(u, v)) obj2 += (IloNum(rna_.get_pij(u, v)) * y(u, v));
                }
            }
            break;
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

void MOIP::define_problem_constraints(string& source)
{

    // ensure there only is 0 or 1 pairing by nucleotide:
    if (verbose_) cout << "\t> ensuring there are at most 1 pairing by nucleotide..." << endl;
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

    // Ensure that the stacking of (i,j) and (i+1,j-1) exists if and only if the pairing (i,j) and (i+1, j-1) exist
    if (verbose_) cout << "\t> ensuring that the stacks are correct..." << endl;
    for (u = 0; u < n - 5; u++) {
        for (v = u + 4; v < n; v++) {
            if (allowed_basepair(u, v) and allowed_basepair(u + 1, v - 1)) {
                IloExpr c7_1(env_);
                IloExpr c7_2(env_);

                c7_1 += y(u, v) + y(u + 1, v - 1);
                c7_2 += y(u, v) + y(u + 1, v - 1) - IloNum(1);

                model_.add(IloNum(2) * x(u, v) <= c7_1);
                if (verbose_) cout << "\t\t" << (2 * x(u,v) <= c7_1) << endl;
                model_.add(x(u, v) >= c7_2);
                if (verbose_) cout << "\t\t" << (x(u, v) >= c7_2) << endl;
            }
        }
    }

    // forbid lonely basepairs if databases other than CaRNAval are being used
    if (source != "rinfolder" and source != "jsonfolder")
    {
        if (verbose_) cout << "\t> forbidding lonely basepairs..." << endl;
        for (u = 0; u < n - 5; u++)
            for (v = u + 4; v < n; v++)
            {
                if (allowed_basepair(u, v))
                {
                    IloExpr c2(env_);
                    c2 += -y(u, v);
                    if (allowed_basepair(u - 1, v + 1)) c2 += y(u - 1, v + 1);
                    if (allowed_basepair(u + 1, v - 1)) c2 += y(u + 1, v - 1);
                    model_.add(c2 >= 0);
                    if (verbose_) cout << "\t\t" << (c2 >= 0) << endl;
                }
            }
    }

    // Forbid pairings inside every motif component if included
    if (verbose_) {
        if (source == "jsonfolder" or source == "rinfolder") {
            cout << "\t> forbidding basepairs inside a motif component if they are not explicitly part of the motif..." << endl;
        } else {
            cout << "\t> forbidding basepairs inside included motif's components..." << endl;
        }
    }
    for (size_t i = 0; i < insertion_sites_.size(); i++)
    {
        Motif& x = insertion_sites_[i];

        for (size_t j = 0; j < x.comp.size(); j++)
        {
            Component& c = x.comp[j];
            IloExpr    c3(env_);
            IloNum     kxi = IloNum(c.k);
            
            uint count = 0;
            if (source == "jsonfolder" or source == "rinfolder") {
                c3 += kxi * C(i, j);
                for (u = c.pos.first ; u <= c.pos.second; u++)
                    for (v = 0; v < n; v++)
                    {
                        if (allowed_basepair(u,v))
                        {
                            bool is_link = false;
                            for (Link link : x.links_)
                                if ((u == link.nts.first and v == link.nts.second) or (u == link.nts.second and v == link.nts.first))
                                {
                                    is_link = true;
                                    break;
                                }

                            if (!is_link)
                            {
                                c3 += y(u, v);
                                count++;
                            }
                        
                        }
                    }
            } else {
                c3 += (kxi - IloNum(2)) * C(i, j);
                for (u = c.pos.first + 1; u < c.pos.second; u++)
                    for (v = 0; v < n; v++)
                    {
                        if (allowed_basepair(u,v))
                        {
                            c3 += y(u, v);
                            count++;
                        }
                    }
            }
            if (count > 0)
            {   
                if (verbose_) cout << "\t\t";
                if (verbose_) cout << x.get_identifier() << '-' << j << ": ";
                if(source == "jsonfolder" or source == "rinfolder") {
                    if (verbose_) cout << (c3 <= kxi) << endl;
                    model_.add(c3 <= kxi);
                } else { 
                    if (verbose_) cout << (c3 <= (kxi - IloNum(2))) << endl;
                    model_.add(c3 <= (kxi - IloNum(2)));
                }
            }
        }
    }
    // Forbid component overlap
    if (verbose_) cout << "\t> forbidding component overlap..." << endl;
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
    if (verbose_) cout << "\t> ensuring that motives cannot be partially included..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++) {
        Motif& x = insertion_sites_[i];
        /*for (uint ii = 0; ii < x.comp.size(); ii++) {
            cout << endl << "insertion (" << i << "): " << x.comp[ii].pos.first << ";" << x.comp[ii].pos.second << endl << endl;
        }*/
        if (x.comp.size() == 1)    // This constraint is for multi-component motives.
            continue;
        IloExpr c5(env_);
        IloNum  jm1 = IloNum(x.comp.size() - 1);
        for (size_t j = 1; j < x.comp.size(); j++) {
            c5 += C(i, j);
        }
        model_.add(c5 == jm1 * C(i, 0));
        if (verbose_) cout << "\t\t> motif " << i << " : " << (c5 == jm1 * C(i, 0)) << endl;
    }

    // basepairs between components
    if (verbose_) cout << "\t> forcing basepairs imposed by a module insertion..." << endl;
    if (source == "jsonfolder" or source == "rinfolder") {
        for (size_t i=0; i < insertion_sites_.size(); i++) {
            Motif&  x   = insertion_sites_[i];
            if (verbose_) cout << "\t\t> motif " << i << " " << x.pos_string() << " (" << x.links_.size() << " canonical pairs)";
            
            for (size_t j=0; j < x.comp.size(); j++) {
                Component& c = x.comp[j];
                IloExpr    c6(env_);
                uint ax(0);

                for (size_t k=0; k < x.links_.size(); k++) // iterate on the motif links
                {
                    size_t ntA = x.links_[k].nts.first;
                    size_t ntB = x.links_[k].nts.second;

                    //check if the component we are in (j) is the first to be linked in the k link
                    if( c.pos.first <= ntA && ntA <= c.pos.second )
                    {
                        if (allowed_basepair(ntA, ntB))
                        {
                            c6 += y(ntA, ntB);
                            ax++;
                        }
                        else // a link is unauthorized, the component cannot be inserted
                        {
                            ax++;
                            break;
                        }       
                    }
                }

                if (ax > 0) {
                    model_.add(c6 >= IloNum(ax) * C(i, j));
                    if (verbose_) cout << "\t\t" << (IloNum(ax) * C(i, j) <= c6);
                }
            }
            if (verbose_) cout << endl;
        }
    }
    else
    {
        // Force basepairs between the end of a component and the beginning of the next
        for (size_t i = 0; i < insertion_sites_.size(); i++)
        {
            Motif&  x   = insertion_sites_[i];
            IloExpr c6p = IloExpr(env_);

            if (allowed_basepair(x.comp[0].pos.first, x.comp.back().pos.second))
                c6p += y(x.comp[0].pos.first, x.comp.back().pos.second);

            if (verbose_) cout << "\t\t" << (IloNum(1) * C(i, 0) <= c6p) << endl;

            model_.add(C(i, 0) <= c6p);

            if (x.comp.size() == 1)    // This constraint is for multi-component motives.
                continue;

            for (size_t j = 0; j < x.comp.size() - 1; j++)
            {
                IloExpr c6 = IloExpr(env_);

                if (allowed_basepair(x.comp[j].pos.second, x.comp[j + 1].pos.first)) //nt u et v
                    c6 += y(x.comp[j].pos.second, x.comp[j + 1].pos.first);

                model_.add(C(i, j) <= c6);

                if (verbose_) cout << "\t\t" << (IloNum(1) * C(i, j) <= c6) << endl;
            }
        }
    }
    
    // Forbid pseudoknots
    if (!this->allow_pk_) {
        if (verbose_) cout << "\t> forbidding pseudoknots..." << endl;
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
        if (verbose_) cout << "\t> Failed to optimize LP: no more solutions to find." << endl;
        // Removing the objective from the model_
        model_.remove(obj);
        model_.remove(bounds);
        return SecondaryStructure(true);
    }
    if (verbose_)
        cout << "\t> Solution status: objective values (" << cplex_.getValue(obj1) << ", " << cplex_.getValue(obj2) << ')';

    // Build a secondary Structure
    SecondaryStructure best_ss = SecondaryStructure(rna_);
    // if (verbose_) cout << "\t\t>retrieveing motifs inserted in the result secondary structure..." << endl;
    for (size_t i = 0; i < insertion_sites_.size(); i++)
        // A constraint requires that all the components are inserted or none, so testing the first is enough:
        if (cplex_.getValue(insertion_dv_[index_of_first_components[i]]) > 0.5) {
            best_ss.insert_motif(insertion_sites_[i]);
        }

    // if (verbose_) cout << "\t\t>retrieving basepairs of the result secondary structure..." << endl;
    for (size_t u = 0; u < rna_.get_RNA_length() - 6; u++)
        for (size_t v = u + 4; v < rna_.get_RNA_length(); v++)
            if (allowed_basepair(u, v))
                if (cplex_.getValue(y(u, v)) > 0.5) {
                    best_ss.set_basepair(u, v);
                }

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

void MOIP::search_between(double lambdaMin, double lambdaMax)
{
    //if (fabs(lambdaMin - lambdaMax) < MOIP::precision_) return;
    if (lambdaMin - lambdaMax > 0.0) return;
    SecondaryStructure s = solve_objective(obj_to_solve_, lambdaMin + MOIP::precision_, lambdaMax);
    //cout << "min: " << lambdaMin << " max: " << lambdaMax << endl;
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
                        cout << "\t> removing structure from Pareto set, obj " << 3 - obj_to_solve_ << " = "
                             << x->get_objective_score(3 - obj_to_solve_) << endl;
                    pareto_.erase(x);
                }
        // search on top
        double min = s.get_objective_score(3 - obj_to_solve_) + precision_;
        double max = lambdaMax;

        if (verbose_)
            cout << std::setprecision(-log10(precision_) + 7) << "\nSolving objective function " << obj_to_solve_
                 << ", on top of " << s.get_objective_score(3 - obj_to_solve_) << ": Obj" << 3 - obj_to_solve_
                 << "  being in [" << std::setprecision(-log10(precision_) + 7) << min << ", "
                 << std::setprecision(-log10(precision_) + 7) << max << "]..." << endl;
        search_between(min, max);


        if (std::abs(max - min) - precision_ > precision_) {

            // search below
            min = lambdaMin;
            max = s.get_objective_score(3 - obj_to_solve_);
            if (verbose_)
                cout << std::setprecision(-log10(precision_) + 7) << "\nSolving objective function " << obj_to_solve_
                     << ", below (or eq. to) " << max << ": Obj" << 3 - obj_to_solve_ << "  being in ["
                     << std::setprecision(-log10(precision_) + 7) << min << ", "
                     << std::setprecision(-log10(precision_) + 7) << max << "]..." << endl;
            search_between(min, max);
        }

    } else {
        if (verbose_) cout << "\t> no solutions found." << endl;
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
    if (verbose_) cout << "\t> adding structure to Pareto set :\t" << s.to_string() << endl; 
    pareto_.push_back(s);
    if (pareto_.size() > max_sol_nbr_) {
        cerr << "\033[31m Quitting because combinatorial issues (>" << max_sol_nbr_ << " solutions in Pareto set). \033[0m" << endl;
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

size_t MOIP::get_xij_index(size_t u, size_t v) const
{
    size_t a, b;
    a = (u < v) ? u : v;
    b = (u > v) ? u : v;
    return index_of_xij_[a][b - 4 - a];
}

void MOIP::remove_solution(uint i) { pareto_.erase(pareto_.begin() + i); }

bool MOIP::allowed_basepair(size_t u, size_t v) const
{
    size_t a, b;
    a = (v > u) ? u : v;
    b = (v > u) ? v : u;
    
    if (b - a < 4) return false;
    if (a >= rna_.get_RNA_length() - 6) return false;
    if (b >= rna_.get_RNA_length()) return false;
    if (get_yuv_index(a, b) == rna_.get_RNA_length() * rna_.get_RNA_length() + 1) {
        //cout << get_yuv_index()
        return false;    // not allowed because proba < theta
    }
    return true;
}

void MOIP::allowed_motifs_from_desc(file_and_mutex arg_struct)
{
    /*
        Searches where to place some DESC module in the RNA
        Too short components are extended in all possible directions.
    */
    path           descfile                 = arg_struct.motif_file;
    mutex&         posInsertionSites_access = arg_struct.posInsertionSites_mutex;

    std::ifstream             motif;
    vector<vector<Component>> vresults;
    string                    line;
    string                    seq;
    vector<string>            component_sequences;
    vector<string>            bases;
    int                       last;
    char                      c    = 'a';
    char*                     prev = &c;
    string                      rna  = rna_.get_seq();

    motif = std::ifstream(descfile.string());
    getline(motif, line);    // ignore "id: number"
    getline(motif, line);    // Bases: 866_G  867_G  868_G  869_G  870_U  871_A ...
    boost::split(bases, line, [prev](char c) {
        bool res = (*prev == ' ' or *prev == ':');
        *prev    = c;
        return (c == ' ' and res);
    });    // get a vector of 866_G, 867_G, etc...

    seq  = "";
    last = stoi(bases[1].substr(0, bases[1].find('_')));
    for (vector<string>::iterator b = bases.begin() + 1; b != bases.end() - 1; b++) {
        char nt  = b->substr(b->find('_') + 1, 1).back();
        int  pos = stoi(b->substr(0, b->find('_')));

        if (pos - last > 5) {    // finish this component and start a new one
            component_sequences.push_back(seq);
            seq = "";
        } else if (pos - last == 2) {
            seq += '.';
        } else if (pos - last == 3) {
            seq += "..";
        } else if (pos - last == 4) {
            seq += "...";
        } else if (pos - last == 5) {
            seq += "....";
        }
        seq += nt;
        last = pos;
    }
    component_sequences.push_back(seq);
    // Now component_sequences is a vector of sequences like {AGCGC, CGU..GUUU}

    // identify components of length 1 or 2, then extend them to length 3
    vector<uint> comp_of_size_1;
    vector<uint> comp_of_size_2;
    for (uint p = 0; p < component_sequences.size(); ++p) {
        if (component_sequences[p].length() == 1) comp_of_size_1.push_back(p);
        if (component_sequences[p].length() == 2) comp_of_size_2.push_back(p);
    }
    if (comp_of_size_1.size() or comp_of_size_2.size()) {
        // We have short components to extend.
        // We will look at all the possible extensions of the motifs for which all 
        // components have length 3 or more, and store the variants in motif_variants:
        vector<vector<string>> motif_variants;

        component_sequences.clear();    // rebuild from scratch
        motif_variants.push_back(component_sequences);
        uint actual_comp = 0;

        seq  = "";
        last = stoi(bases[1].substr(0, bases[1].find('_')));
        for (vector<string>::iterator b = bases.begin() + 1; b < bases.end() - 1; b++) {
            int  pos = stoi(b->substr(0, b->find('_')));
            char nt  = b->substr(b->find('_') + 1, 1).back();
            if (comp_of_size_1.size() and actual_comp == comp_of_size_1[0])    // we are on the first component of size 1
            {
                b--;
                nt          = b->substr(b->find('_') + 1, 1).back();
                string seq1 = "";
                seq1 += nt;
                seq1 += "..";
                string seq2 = ".";
                seq2 += nt;
                seq2 += ".";
                string seq3 = "..";
                seq3 += nt;
                uint end = motif_variants.size();    // before to add the new ones
                for (uint u = 0; u < end; ++u) {
                    motif_variants.push_back(motif_variants[u]);    // copy 1 for seq2
                    motif_variants.back().push_back(seq2);
                    motif_variants.push_back(motif_variants[u]);    // copy 2 for seq3
                    motif_variants.back().push_back(seq3);
                    motif_variants[u].push_back(seq1);
                }
                seq = "";
                actual_comp++;
                comp_of_size_1.erase(comp_of_size_1.begin());    // the first element has been processed, remove it
                last = pos;
            } else if (comp_of_size_2.size() and actual_comp == comp_of_size_2[0]) {    // we are on the first component of size 2
                b--;
                nt = b->substr(b->find('_') + 1, 1).back();
                b++;    // skip the next nucleotide
                char next   = b->substr(b->find('_') + 1, 1).back();
                last        = stoi(b->substr(0, b->find('_')));
                string seq1 = "";
                seq1 += nt;
                seq1 += next;
                seq1 += ".";
                string seq2 = ".";
                seq2 += nt;
                seq2 += next;
                uint end = motif_variants.size();    // before to add the new one
                for (uint u = 0; u < end; ++u) {
                    motif_variants.push_back(motif_variants[u]);    // copy 1 for seq2
                    motif_variants.back().push_back(seq2);
                    motif_variants[u].push_back(seq1);
                }
                seq = "";
                actual_comp++;
                comp_of_size_2.erase(comp_of_size_2.begin());    // the first element has been processed, remove it
            } else {                                             // we are on a longer component
                if (pos - last > 5) {                            // finish this component and start a new one
                    actual_comp++;
                    for (vector<string>& c_s : motif_variants) c_s.push_back(seq);
                    seq = "";
                } else if (pos - last == 2) {
                    seq += '.';
                } else if (pos - last == 3) {
                    seq += "..";
                } else if (pos - last == 4) {
                    seq += "...";
                } else if (pos - last == 5) {
                    seq += "....";
                }
                seq += nt;
                last = pos;
            }
        }
        for (auto c_s : motif_variants)
            if (seq.length()) c_s.push_back(seq);    // pushing the last one after iterating over the bases

        // We need to search for the different positions where to insert the first component
        for (auto c_s : motif_variants) {
            vector<vector<Component>> new_results = find_next_ones_in(rna, 0, c_s);
            vresults.insert(vresults.end(), new_results.begin(), new_results.end());
        }

    } 
    else 
    {
        // No multiple motif variants : we search in a single vector component_sequences
        // We need to search for the different positions where to insert the first component
        vresults = find_next_ones_in(rna, 0, component_sequences);
    }

    // Now create proper motifs with Motif class
    for (vector<Component>& v : vresults) {
        Motif temp_motif = Motif(v, path(descfile).stem().string());

        // Check if the probabilities allow to keep this Motif:
        bool unprobable = false;
        if (!allowed_basepair(temp_motif.comp[0].pos.first, temp_motif.comp.back().pos.second))
            unprobable = true;
        for (size_t j = 0; j < temp_motif.comp.size() -1; j++)
            if (!allowed_basepair(temp_motif.comp[j].pos.second, temp_motif.comp[j+1].pos.first))
                unprobable = true;
        if (unprobable) continue;

        // Add it to the results vector
        unique_lock<mutex> lock(posInsertionSites_access);
        insertion_sites_.push_back(temp_motif);
        lock.unlock();
    }
}

void MOIP::allowed_motifs_from_rin(file_and_mutex arg_struct)
{
    // Searches where to place some RINs in the RNA

    path           rinfile                  = arg_struct.motif_file;
    mutex&         posInsertionSites_access = arg_struct.posInsertionSites_mutex;

    std::ifstream                 motif;
    string                           filepath = rinfile.string();
    vector<vector<Component>>     vresults, r_vresults;
    vector<string>                component_sequences;
    uint                          carnaval_id;
    string                        line, filenumber;
    string                        rna = rna_.get_seq();
    string                        reversed_rna = rna_.get_seq();

    std::reverse(reversed_rna.begin(), reversed_rna.end());
    filenumber = filepath.substr(filepath.find("Subfiles/")+9, filepath.find(".txt"));
    carnaval_id = 1 + stoi(filenumber); // Start counting at 1 to be consistant with the website numbering

    motif = std::ifstream(rinfile.string());
    getline(motif, line); //skip the header_link line
    getline(motif, line); //get the links line
    getline(motif, line); //skip the header_comp line
    // std::cout << "RIN " << carnaval_id << ": \t";
    while (getline(motif, line))
    {
        // lines are formatted like:
        // pos;k;seq
        // 0,1;2;GU
        if (line == "\n") break; //skip last line (empty)
        size_t index = line.find(';', line.find(';') + 1); // find the second ';'
        component_sequences.push_back(line.substr(index+1, string::npos)); // new component sequence
        // std::cout << line.substr(index+1, string::npos) << " ";
    }

    vresults     = find_next_ones_in(rna, 0, component_sequences);

    for (vector<Component>& v : vresults)
    {
        // Now create an actual Motif with the class
        Motif temp_motif = Motif(v, rinfile, carnaval_id, false);

        bool unprobable = false;
        for (const Link& l : temp_motif.links_)
        {
            if (!allowed_basepair(l.nts.first,l.nts.second))
                unprobable = true;
        }
        if (unprobable) continue;

        // Add it to the results vector
        unique_lock<mutex> lock(posInsertionSites_access);
        insertion_sites_.push_back(temp_motif);
        lock.unlock();
    }
}

void MOIP::allowed_motifs_from_json(motif_and_mutex arg_struct)
{
    // Searches where to place some JSON motifs in the RNA

    mutex&                    posInsertionSites_access = arg_struct.posInsertionSites_mutex;
    json_elem                 it = arg_struct.motif;
    vector<vector<Component>> vresults;
    vector<string>            component_sequences;
    string                    subseq, seq, fullseq, struct2d;
    size_t                    end;

    // Retrieve sequenc eand structure
    seq = it.value()["sequence"];
    fullseq = string(seq);
    struct2d = it.value()["struct2d"];
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

    // Now split the sequence into components
    while(seq.find("&") != string::npos) {
        end = seq.find("&");
        subseq = seq.substr(0, end);
        seq = seq.substr(end + 1);
        component_sequences.push_back(subseq); // new component sequence
    } 
    if (!seq.empty()) {
        component_sequences.push_back(seq);
    }

    // Place the components by pattern-matching
    vresults = find_next_ones_in(rna_.get_seq(), 0, component_sequences);
    
    // Check if the results are possible (likeliness)
    for (vector<Component>& v : vresults)
    {
        
        // Now create proper Motifs with Motif class
        Motif temp_motif = Motif(v, it.key(), struct2d);

        if (verbose_) cout << "\t> Considering motif JSON " << temp_motif.pos_string() << "\t" 
                           << fullseq << ", " << struct2d << " ";

        // Check if the motif can be inserted, checking the basepairs probabilities and theta
        bool unprobable = false;
        if (!temp_motif.links_.size()) {
            if (verbose_) cout << "discarded, no constraints on the secondary structure, it is a useless motif." << endl;
            continue;
        }
        if (verbose_) cout << " + links at positions ";
        for (const Link& l : temp_motif.links_)
        {
            if (verbose_) cout << l.nts.first << ',' << l.nts.second << ' '; 
            if (!allowed_basepair(l.nts.first, l.nts.second)) {
                if (verbose_) cout << "(unlikely) ";
                unprobable = true;
            }
        }
        if (unprobable) {
            if (verbose_) cout << ", discarded because of unlikely or impossible basepairs" << endl;
            continue;
        }
        if (verbose_) cout << endl;

        // Add it to the results vector
        unique_lock<mutex> lock(posInsertionSites_access);
        insertion_sites_.push_back(temp_motif);
        lock.unlock();
    }
}