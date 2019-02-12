/***
        Biominserter, Louis Becquey, nov 2018
***/

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

#include "MOIP.h"
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

Motif parse_csv_line(string line)
{
    vector<string> tokens;
    boost::split(tokens, line, boost::is_any_of(","));
    Motif m;
    m.atlas_id = tokens[0];
    m.score    = stoi(tokens[2]);
    m.comp.push_back(Component(make_pair<int, int>(stoi(tokens[3]), stoi(tokens[4]))));
    if (tokens[5] != "-") m.comp.push_back(Component(make_pair<int, int>(stoi(tokens[5]), stoi(tokens[6]))));
    m.reversed = (tokens[1] == "True");
    return m;
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
    string             line;
    ifstream           motifs;
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
    motifs = ifstream(csvname);
    getline(motifs, line);    // skip header
    while (getline(motifs, line)) posInsertionSites.push_back(parse_csv_line(line));
    if (verbose)
        cout << "\t>" << csvname << " successfuly loaded (" << posInsertionSites.size() << " insertion sites)" << endl;

    /*  FIND K-PARETO SETS  */

    MOIP myMOIP = MOIP(myRNA, posInsertionSites, 1, theta_p_threshold, verbose);
    try {
        bestSSO1 = myMOIP.solve_objective(1, -__DBL_MAX__, __DBL_MAX__);
        bestSSO2 = myMOIP.solve_objective(2, -__DBL_MAX__, __DBL_MAX__);
        bestSSO1.set_pareto_set(1);
        bestSSO2.set_pareto_set(1);
        myMOIP.add_solution(bestSSO1);
        myMOIP.add_solution(bestSSO2);
        if (verbose) {
            cout << endl << "Best solution according to objective 1 :" << bestSSO1.to_string() << endl;
            cout << "Best solution according to objective 2 :" << bestSSO2.to_string() << endl;
        }

        // extend to the whole pareto set
        if (MOIP::obj_to_solve_ == 1) {
            if (verbose) cout << endl << "Solving obj1 on top of best solution 1." << endl;
            myMOIP.search_between(bestSSO1.get_objective_score(2) + MOIP::epsilon_, bestSSO2.get_objective_score(2));
            if (verbose) cout << endl << "Solving obj1 below best solution 1." << endl;
            myMOIP.search_between(-__DBL_MAX__, bestSSO1.get_objective_score(2));
        } else {
            if (verbose) cout << endl << "Solving obj2 on top of best solution 2." << endl;
            myMOIP.search_between(bestSSO2.get_objective_score(1) + MOIP::epsilon_, bestSSO1.get_objective_score(1));
            if (verbose) cout << endl << "Solving obj2 below best solution 2." << endl;
            myMOIP.search_between(-__DBL_MAX__, bestSSO2.get_objective_score(1));
        }
    } catch (IloAlgorithm::NotExtractedException& e) {
        cerr << e << endl;
        exit(EXIT_FAILURE);
    } catch (IloCplex::Exception& e) {
        cerr << e << endl;
        exit(EXIT_FAILURE);
    }

    /*  REMOVE SOLUTIONS WITH TOO HIGH LABEL  */

    vector<size_t> to_remove;
    if (verbose) cout << endl;
    for (uint i = 0; i < myMOIP.get_n_solutions(); i++)
        if (myMOIP.solution(i).get_pareto_set() > 1) {    // Some solution is fromm a Pareto set of too high order
            if (verbose)
                cout << "Removing structure from Pareto set " << myMOIP.solution(i).get_pareto_set() << " : "
                     << myMOIP.solution(i).to_string() << endl;
            to_remove.push_back(i);
        }
    if (to_remove.size()) {
        for (size_t i = to_remove.size() - 1; i != 0; i--) myMOIP.remove_solution(to_remove[i]);
        myMOIP.remove_solution(to_remove[0]);
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
