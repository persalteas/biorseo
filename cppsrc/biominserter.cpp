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
    if (argc != 5) {
        cerr << argc << " arguments specified !" << endl;
        cerr << "Please specify the following input files:" << endl;
        cerr << "biominserter sequence.fasta insertion.sites.csv prob_threshold verbose" << endl;
        return EXIT_FAILURE;
    }


    const char* inputName         = argv[1];
    const char* csvname           = argv[2];
    bool        verbose           = (atoi(argv[4]) != 0);
    string      fastaname         = string(inputName);
    string      basename          = remove_ext(inputName, '.', '/');
    float       theta_p_threshold = atof(argv[3]);
    if (verbose) cout << "Reading input files..." << endl;
    if (access(inputName, F_OK) == -1) {
        cerr << inputName << " not found" << endl;
        return EXIT_FAILURE;
    }
    if (access(csvname, F_OK) == -1) {
        cerr << csvname << " not found" << endl;
        return EXIT_FAILURE;
    }

    // load fasta file
    list<Fasta> f;
    Fasta::load(f, inputName);
    list<Fasta>::iterator fa = f.begin();
    if (verbose) cout << "loading " << fa->name() << "..." << endl;
    RNA myRNA = RNA(fa->name(), fa->seq(), verbose);
    if (verbose) cout << "\t>" << inputName << " successfuly loaded (" << myRNA.get_RNA_length() << " nt)" << endl;

    // load CSV file
    string   line;
    ifstream motifs = ifstream(csvname);
    getline(motifs, line);    // skip header
    vector<Motif> posInsertionSites;
    while (getline(motifs, line)) {
        posInsertionSites.push_back(parse_csv_line(line));
    }
    if (verbose)
        cout << "\t>" << csvname << " successfuly loaded (" << posInsertionSites.size() << " insertion sites)" << endl;

    // creating the Multi-Objective problem:
    MOIP myMOIP = MOIP(myRNA, posInsertionSites, theta_p_threshold, verbose);
    // finding the best SecondaryStructures for each objective
    try {
        SecondaryStructure bestSSO1 = myMOIP.solve_objective(1, -__DBL_MAX__, __DBL_MAX__);
        SecondaryStructure bestSSO2 = myMOIP.solve_objective(2, -__DBL_MAX__, __DBL_MAX__);
        if (verbose) cout << "Best solution according to objective 1 :" << bestSSO1.to_string() << endl;
        if (verbose) cout << "Best solution according to objective 2 :" << bestSSO2.to_string() << endl;

        // extend to the whole pareto set
        myMOIP.add_solution(bestSSO1);
        myMOIP.add_solution(bestSSO2);
        myMOIP.extend_pareto(bestSSO2.get_objective_score(1), bestSSO1.get_objective_score(1));

        // print the pareto set
        if (verbose) {
            cout << endl << endl << "---------------------------------------------------------------" << endl;
            cout << "Whole Pareto Set:" << endl;
            for (uint i = 0; i < myMOIP.get_n_solutions(); i++) {
                myMOIP.solution(i).print();
            }
            cout << endl;

            cout << posInsertionSites.size() << " candidate insertion sites, " << myMOIP.get_n_solutions()
                 << " solutions kept." << endl;
            cout << "Best value for Motif insertion objective: " << bestSSO1.get_objective_score(1) << endl;
            cout << "Best value for structure expected accuracy: " << bestSSO2.get_objective_score(2) << endl;
        }
    } catch (IloAlgorithm::NotExtractedException& e) {
        cerr << e << endl;
        exit(EXIT_FAILURE);
    } catch (IloCplex::Exception& e) {
        cerr << e << endl;
        exit(EXIT_FAILURE);
    }

    // Save the pareto set to a file
    if (verbose) cout << "Saving structures to " << basename << ".biom..." << endl;
    ofstream outfile;
    outfile.open(basename + ".biom");
    outfile << fa->name() << endl << fa->seq() << endl;
    for (uint i = 0; i < myMOIP.get_n_solutions(); i++) {
        outfile << myMOIP.solution(i).to_string() << endl;
    }
    outfile.close();

    return EXIT_SUCCESS;
}
