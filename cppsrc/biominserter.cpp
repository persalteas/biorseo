/***
        Biominserter, Louis Becquey, nov 2018
***/

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

#include "MOIP.h"
#include "Motif.h"
#include "fa.h"

using namespace std;

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
    /*  ARGUMENT CHECKING  */

    if (argc != 6) {
        cerr << argc << " arguments specified !" << endl;
        cerr << "Please specify the following input files:" << endl;
        cerr << "biominserter sequence.fasta insertion.sites.csv prob_threshold verbose obj" << endl;
        return EXIT_FAILURE;
    }

    /*  VARIABLE DECLARATIONS  */

    const char*        inputName         = argv[1];
    const char*        csvname           = argv[2];
    bool               verbose           = (atoi(argv[4]) != 0);
    string             basename          = remove_ext(inputName, '.', '/');
    float              theta_p_threshold = atof(argv[3]);
    list<Fasta>        f;
    vector<Motif>      posInsertionSites;
    ofstream           outfile;
    SecondaryStructure bestSSO1, bestSSO2;
    RNA                myRNA;
    MOIP::obj_to_solve_ = atoi(argv[5]);

    /*  FILE PARSING  */

    // load fasta file
    if (verbose) cout << "Reading input files..." << endl;
    if (access(inputName, F_OK) == -1) {
        cerr << inputName << " not found" << endl;
        return EXIT_FAILURE;
    }
    Fasta::load(f, inputName);
    list<Fasta>::iterator fa = f.begin();
    if (verbose) cout << "loading " << fa->name() << "..." << endl;
    myRNA = RNA(fa->name(), fa->seq(), verbose);
    if (verbose) cout << "\t>" << inputName << " successfuly loaded (" << myRNA.get_RNA_length() << " nt)" << endl;

    // load CSV file
    if (access(csvname, F_OK) == -1) {
        cerr << csvname << " not found" << endl;
        return EXIT_FAILURE;
    }
    posInsertionSites = load_desc_folder(csvname);
    if (verbose)
        cout << "\t>" << csvname << " successfuly loaded (" << posInsertionSites.size() << " insertion sites)" << endl;

    /*  FIND PARETO SET  */

    MOIP               myMOIP = MOIP(myRNA, posInsertionSites, theta_p_threshold, verbose);
    double             min, max;
    IloConstraintArray F(myMOIP.get_env());

    try {
        bestSSO1 = myMOIP.solve_objective(1, -__DBL_MAX__, __DBL_MAX__);
        bestSSO2 = myMOIP.solve_objective(2, -__DBL_MAX__, __DBL_MAX__);
        if (verbose) {
            cout << endl << "Best solution according to objective 1 :" << bestSSO1.to_string() << endl;
            cout << "Best solution according to objective 2 :" << bestSSO2.to_string() << endl;
        }

        // extend the Pareto set on top
        if (MOIP::obj_to_solve_ == 1) {
            myMOIP.add_solution(bestSSO1);
            min = bestSSO1.get_objective_score(2) + MOIP::precision_;
            max = bestSSO2.get_objective_score(2);
            if (verbose) cout << endl << "Solving obj1 on top of best solution 1." << endl;
        } else {
            myMOIP.add_solution(bestSSO2);
            min = bestSSO2.get_objective_score(1) + MOIP::precision_;
            max = bestSSO1.get_objective_score(1);
            if (verbose) cout << endl << "Solving obj2 on top of best solution 2." << endl;
        }

        if (verbose)
            cout << std::setprecision(-log10(MOIP::precision_) + 4) << "\nSolving objective function "
                 << MOIP::obj_to_solve_ << ", on top of " << min << ": Obj" << 3 - MOIP::obj_to_solve_ << "  being in ["
                 << min << ", " << max << "]..." << endl;
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
            cout << std::setprecision(-log10(MOIP::precision_) + 4) << "\nSolving objective function "
                 << MOIP::obj_to_solve_ << ", below (or eq. to) " << max << ": Obj" << 3 - MOIP::obj_to_solve_
                 << "  being in [" << min << ", " << max << "]..." << endl
                 << "\t>forbidding " << F.getSize() << " solutions found in [" << std::setprecision(10)
                 << min - MOIP::precision_ << ", " << max + MOIP::precision_ << ']' << endl;
        myMOIP.search_between(min, max);

    } catch (IloCplex::Exception& e) {
        cerr << "Cplex Exception: " << e.getMessage() << endl;
        exit(EXIT_FAILURE);
    }

    /*  DISPLAY RESULTS  */

    // print the pareto set
    if (verbose) {
        cout << endl << endl << "---------------------------------------------------------------" << endl;
        cout << "Whole Pareto Set:" << endl;
        for (uint i = 0; i < myMOIP.get_n_solutions(); i++) myMOIP.solution(i).print();
        cout << endl;
        cout << posInsertionSites.size() << " candidate insertion sites, " << myMOIP.get_n_solutions() << " solutions kept." << endl;
        cout << "Best value for Motif insertion objective: " << bestSSO1.get_objective_score(1) << endl;
        cout << "Best value for structure expected accuracy: " << bestSSO2.get_objective_score(2) << endl;
    }

    // Save it to file
    if (verbose) cout << "Saving structures to " << basename << ".biom..." << endl;
    outfile.open(basename + ".biom");
    outfile << fa->name() << endl << fa->seq() << endl;
    for (uint i = 0; i < myMOIP.get_n_solutions(); i++) outfile << myMOIP.solution(i).to_string() << endl;
    outfile.close();

    /*  QUIT  */

    return EXIT_SUCCESS;
}
