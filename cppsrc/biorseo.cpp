/***
        Biorseo, Louis Becquey, nov 2018-nov 2021
        Special thanks to Lénaic Durand and Nathalie Bertrand for working a lot on version 2.0
        louis.becquey@univ-evry.fr (likely expired email account, find me on the internet)
***/
#include <cmath> 

#include <algorithm>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "MOIP.h"
#include "Motif.h"
#include "fa.h"

using namespace std;
namespace po = boost::program_options;

string remove_ext(const char* mystr, char dot, char sep)
{
    // COPYPASTA from stackoverflow

    char *retstr, *lastdot, *lastsep;

    // Error checks and allocate string.
    if (mystr == nullptr) return nullptr;
    if ((retstr = static_cast<char*>(malloc(strlen(mystr) + 1))) == nullptr) return nullptr;

    // Make a copy and find the relevant characters.
    strcpy(retstr, mystr);
    lastdot = strrchr(retstr, dot);
    lastsep = (sep == 0) ? nullptr : strrchr(retstr, sep);

    // If it has an extension separator.
    if (lastdot != nullptr) {
        // and it's before the extenstion separator.
        if (lastsep != nullptr) {
            if (lastsep < lastdot) {
                // then remove it.
                *lastdot = '\0';
            }
        } else {
            // Has extension separator with no path separator.
            *lastdot = '\0';
        }
    }

    // Return the modified string.
    return string(retstr);
}

int main(int argc, char* argv[])
{
    /*  VARIABLE DECLARATIONS  */

    string             inputName, outputName, motifs_path_name, basename;
    bool               verbose = false;
    float              theta_p_threshold;
    char               obj_function_nbr;
    char               mea_or_mfe = 'b'; // a for MFE, b for MEA
    list<Fasta>        f;
    ofstream           outfile;
    SecondaryStructure bestSSO1, bestSSO2;
    RNA                myRNA;

    /*  ARGUMENT CHECKING  */
    
    po::options_description desc("Options");
    desc.add_options()
    ("help,h", "Print the help message")
    ("version", "Print the program version")
    ("seq,s", po::value<string>(&inputName)->required(), "Fasta file containing the RNA sequence")
    ("descfolder,d", po::value<string>(&motifs_path_name), "A folder containing modules in .desc format, as produced by Djelloul & Denise's catalog program (deprecated)")
    ("rinfolder,r", po::value<string>(&motifs_path_name), "A folder containing CaRNAval's RINs in .txt format, as produced by script transform_caRNAval_pickle.py")
    ("jsonfolder,j", po::value<string>(&motifs_path_name), "A folder containing a custom motif library in .json format")
    ("pre-placed,x", po::value<string>(&motifs_path_name), "A CSV file providing motif insertion sites obtained with another tool.")
    ("function,f", po::value<char>(&obj_function_nbr)->default_value('B'), 
                    "(A, B, C, D, E or F) Objective function to score module insertions:\n"
                    "  (A) insert big modules\n  (B) light, high-order modules\n"
                    "  (C) well-scored modules\n  (D) light, high-order, well-scored\n    modules\n"
                    "  (E, F) insert big modules with many\n   contacts with proteins, different\n   ponderations.\n"
                    "  C and D require position scores\n  provided by --pre-placed.\n"
                    "  E and F require protein-contact\n  information and should be\n  used only with --jsonfolder.")    
    ("mfe,E", "Minimize stacking energies\n  (leads to MFE extimator)")
    ("mea,A", "(default) Maximize expected accuracy\n  (leads to MEA estimator)")
    ("first-objective,c",   po::value<unsigned int>(&MOIP::obj_to_solve_)->default_value(2), 
                    "(1 or 2) Objective to solve in the mono-objective portions of the algorithm.\n"
                    "  (1) is the module objective,\n  given by --function\n"
                    "  (2) is energy-based objective,\n  either MFE or MEA")
    ("output,o",            po::value<string>(&outputName), "A file to summarize the computation results")
    ("theta,t",             po::value<float>(&theta_p_threshold)->default_value(1e-3, "0.001"), "Pairing probability threshold to consider or not the possibility of pairing")
    ("disable-pseudoknots,n", "Add constraints forbidding the formation of pseudoknots")
    ("limit,l", po::value<unsigned int>(&MOIP::max_sol_nbr_)->default_value(500), "Intermediate number of solutions in the Pareto set above which we give up the calculation.")
    ("verbose,v", "Print what is happening to console");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    basename = remove_ext(inputName.c_str(), '.', '/');

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);    // can throw

        if (vm.count("help") or vm.count("-h")) {
            cout << "Biorseo, Bi-Objective RNA Structure Efficient Optimizer" << endl
                << "Bio-objective integer linear programming framework to predict RNA secondary structures by including known RNA modules." << endl
                << "Developped by Louis Becquey, 2018-2021\nLénaïc Durand, 2019\nNathalie Bernard, 2021" << endl << endl
                << "Usage:\tYou must provide:\n\t1) a FASTA input file with -s," << endl
                << "\t2) a module type with --rna3dmotifs, --carnaval, --contacts or --pre-placed," << endl
                << "\t3) one module-based scoring function with --func A, B, C, D, E or F," << endl
                << "\t4) one energy-based scoring function with --mfe or --mea," << endl
                << "\t5) how to display results: in console (-v), or in a result file (-o)." << endl
                << endl
                << desc << endl;
            return EXIT_SUCCESS;
        }
        if (vm.count("version")) {
            cout << "Biorseo v2.1, dockerized, November 2021" << endl;
            return EXIT_SUCCESS;
        }
        po::notify(vm);    // throws on error, so do after help in case there are any problems
        if (vm.count("mfe")) mea_or_mfe = 'a';
        if (vm.count("mea")) mea_or_mfe = 'b';
        if (vm.count("verbose")) verbose = true;
        if (vm.count("disable-pseudoknots")) MOIP::allow_pk_ = false;
    } catch (po::error& e) {
        cerr << "ERROR: \033[31m" << e.what() << "\033[0m" << endl;
        cerr << desc << endl;
        return EXIT_FAILURE;
    }

    MOIP::obj_function_nbr_ = obj_function_nbr;
    MOIP::obj_function2_nbr_ = mea_or_mfe;

    /*  FILE PARSING  */

    // load fasta file
    if (verbose) cout << "Reading input files..." << endl;
    if (access(inputName.c_str(), F_OK) == -1) {
        cerr << "\033[31m" << inputName << " not found\033[0m" << endl;
        return EXIT_FAILURE;
    }
    Fasta::load(f, inputName.c_str());
    list<Fasta>::iterator fa = f.begin();
    if (verbose) cout << "loading " << fa->name() << "..." << endl;
    myRNA = RNA(fa->name(), fa->seq(), verbose);
    if (verbose) cout << "\t> " << inputName << " successfuly loaded (" << myRNA.get_RNA_length() << " nt)" << endl;

    // check motif folder exists
    if (access(motifs_path_name.c_str(), F_OK) == -1) {
        cerr << "\033[31m" << motifs_path_name << " not found\033[0m" << endl;
        return EXIT_FAILURE;
    }


    /*  FIND PARETO SET  */
    string source;
    if (vm.count("rinfolder"))
        source = "rinfolder";
    else if (vm.count("descfolder"))
        source = "descfolder";
    else if (vm.count("jsonfolder"))
        source = "jsonfolder";
    else if (vm.count("pre-placed"))
        source = "csvfile";
    else
        cerr << "ERR: no source of modules provided !" << endl;
    
    MOIP               myMOIP = MOIP(myRNA, source, motifs_path_name.c_str(), theta_p_threshold, verbose);
    double             min, max;
    IloConstraintArray F(myMOIP.get_env());

    if (verbose)
        cout << "Solving..." << endl;
    try {
        bestSSO1 = myMOIP.solve_objective(1, -__DBL_MAX__, __DBL_MAX__);
        if (verbose) cout << endl;
        bestSSO2 = myMOIP.solve_objective(2, -__DBL_MAX__, __DBL_MAX__);
        if (verbose) {
            cout << endl << "Best solution according to objective 1 :" << bestSSO1.to_string() << endl;
            cout << "Best solution according to objective 2 :" << bestSSO2.to_string() << endl;
        }

        // extend the Pareto set on top
        if (MOIP::obj_to_solve_ == 1) {
            myMOIP.add_solution(bestSSO1);
            min = bestSSO1.get_objective_score(2) + MOIP::precision_;
            max = bestSSO2.get_objective_score(2) + MOIP::precision_;
            if (verbose) cout << endl << "Solving obj1 on top of best solution 1." << endl;
        } else {
            myMOIP.add_solution(bestSSO2);
            min = bestSSO2.get_objective_score(1) + MOIP::precision_;
            max = bestSSO1.get_objective_score(1) + MOIP::precision_;
            if (verbose) cout << endl << "Solving obj2 on top of best solution 2." << endl;
        }

        if (verbose)
            cout << setprecision(-log10(MOIP::precision_) + 4) << "\nSolving objective function " << MOIP::obj_to_solve_ << ", on top of "
                 << min << ": Obj" << 3 - MOIP::obj_to_solve_ << "  being in [" << min << ", " << max << "]..." << endl;
        myMOIP.search_between(min, max);

        // extend the Pareto set below
        if (MOIP::obj_to_solve_ == 1) {
            if (verbose) cout << endl << "Solving obj1 below best solution 1." << endl;
            min = -__DBL_MAX__;
            max = bestSSO1.get_objective_score(2);
        } else {
            if (verbose) cout << endl << "Solving obj2 below best solution 2." << endl;
            min = -__DBL_MAX__;
            max = bestSSO2.get_objective_score(1);
        }
        if (verbose)
            cout << setprecision(-log10(MOIP::precision_) + 4) << "\nSolving objective function " << MOIP::obj_to_solve_
                 << ", below (or eq. to) " << max << ": Obj" << 3 - MOIP::obj_to_solve_ << "  being in [" << min << ", "
                 << max << "]..." << endl
                 << "\t> Forbidding " << F.getSize() << " solutions found in [" << setprecision(10)
                 << min - MOIP::precision_ << ", " << max + MOIP::precision_ << ']' << endl;
        myMOIP.search_between(min, max);

    } catch (IloCplex::Exception& e) {
        cerr << "\033[31mCplex Exception: " << e.getMessage() << "\033[0m" << endl;
        exit(EXIT_FAILURE);
    }

    /*  DISPLAY RESULTS  */

    // print the pareto set
    if (verbose) {
        cout << endl << endl << "---------------------------------------------------------------" << endl;
        cout  << "Whole Pareto Set:" << endl;
        for (uint i = 0; i < myMOIP.get_n_solutions(); i++) myMOIP.solution(i).print();
        cout << endl;
        cout << myMOIP.get_n_candidates() << " candidate insertion sites, " << myMOIP.get_n_solutions() << " solutions kept." << endl;
        cout << "Best value for Motif insertion objective: " << bestSSO1.get_objective_score(1) << endl;
        cout << "Best value for structure expected accuracy: " << bestSSO2.get_objective_score(2) << endl;
    }

    // Save it to file
    if (vm.count("output")) {
        if (verbose) cout << "Saving structures to " << outputName << "..." << endl;
        outfile.open(outputName);
        outfile << fa->name() << endl << fa->seq() << endl;
    
        for (uint i = 0; i < myMOIP.get_n_solutions(); i++) {
            outfile << myMOIP.solution(i).to_string() << endl << structure_with_contacts(myMOIP.solution(i)) << endl;
            string str1 = myMOIP.solution(i).to_string();

        }
        outfile.close();
    }

    /*  QUIT  */

    return EXIT_SUCCESS;
}
