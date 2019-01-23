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

Motif parse_csv_line(string line)
{
    vector<string> tokens;
    boost::split(tokens, line, boost::is_any_of(","));
    Motif m;
    m.atlas_id = tokens[0];
    m.comp.push_back(Component(make_pair<int, int>(stoi(tokens[3]), stoi(tokens[4])), stoi(tokens[2])));
    if (tokens[5] != "-")
        m.comp.push_back(Component(make_pair<int, int>(stoi(tokens[5]), stoi(tokens[6])), stoi(tokens[2])));
    m.reversed = (tokens[1] == "True");
    return m;
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        cerr << argc << " arguments specified !" << endl;
        cerr << "Please specify the following input files:" << endl;
        cerr << "biominserter sequence.fasta insertion.sites.csv prob_threshold" << endl;
        return EXIT_FAILURE;
    }


    const char* inputName         = argv[1];
    const char* csvname           = argv[2];
    float       theta_p_threshold = atof(argv[3]);
    cout << "Reading input files..." << endl;
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
    cout << "loading " << fa->name() << "..." << endl;
    RNA myRNA = RNA(fa->name(), fa->seq());
    cout << "\t>" << inputName << " successfuly loaded (" << myRNA.get_RNA_length() << " nt)" << endl;

    // load CSV file
    string   line;
    ifstream motifs = ifstream(csvname);
    getline(motifs, line);    // skip header
    vector<Motif> posInsertionSites;
    while (getline(motifs, line)) {
        posInsertionSites.push_back(parse_csv_line(line));
    }
    cout << "\t>" << csvname << " successfuly loaded (" << posInsertionSites.size() << " insertion sites)" << endl;

    // creating the Multi-Objective problem:
    MOIP myMOIP = MOIP(myRNA, posInsertionSites, theta_p_threshold);
    // finding the best SecondaryStructures for each objective
    try {
        SecondaryStructure bestSSO1 = myMOIP.solve_objective(1, -__DBL_MAX__, __DBL_MAX__);
        SecondaryStructure bestSSO2 = myMOIP.solve_objective(2, -__DBL_MAX__, __DBL_MAX__);
        cout << "Best solution according to objective 1 :" << bestSSO1.to_string() << endl;
        cout << "Best solution according to objective 2 :" << bestSSO2.to_string() << endl;

        // extend to the whole pareto set
        myMOIP.add_solution(bestSSO1);
        myMOIP.extend_pareto(bestSSO1.get_objective_score(2), bestSSO2.get_objective_score(2));

        // print the pareto set
        cout << endl << endl << "---------------------------------------------------------------" << endl;
        cout << "Whole Pareto Set:" << endl;
        for (uint i = 0; i < myMOIP.get_n_solutions(); i++) {
            myMOIP.solution(i).print();
        }
        cout << endl;
        cout << posInsertionSites.size() << " candidate insertion sites, " << myMOIP.get_n_solutions() << " solutions kept." << endl;
        cout << "Best value for Motif insertion objective: " << bestSSO1.get_objective_score(1) << endl;
        cout << "Best value for structure expected accuracy: " << bestSSO2.get_objective_score(2) << endl;
    } catch (IloAlgorithm::NotExtractedException& e) {
        cerr << e << endl;
        exit(EXIT_FAILURE);
    } catch (IloCplex::Exception& e) {
        cerr << e << endl;
        exit(EXIT_FAILURE);
    }


    return EXIT_SUCCESS;
}
