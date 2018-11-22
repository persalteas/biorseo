/***
        Biominserter
        Louis Becquey, starting from Audrey Legendre's code
        contains adapted pieces of code from IPknot and Nupack
        nov 2018
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
#include "nupack.h"

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
    // float   time;
    // clock_t t1, t2;
    // t1 = clock();

    MOIP::ncores = thread::hardware_concurrency() - 1;

    if (argc != 3) {
        cerr << argc << " arguments specified !" << endl;
        cerr << "Please specify the following input files:" << endl;
        cerr << "biominserter sequence.fasta insertion.sites.csv" << endl;
        return EXIT_FAILURE;
    }


    const char* inputName = argv[1];
    const char* csvname   = argv[2];
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
    MOIP myMOIP = MOIP(myRNA, posInsertionSites);    // using the constructor with arguments automatically defines the decision variables.

    // finding the best SecondaryStructures for each objective
    double max = myRNA.get_RNA_length();
    try {
        myMOIP.solve_objective(1, -max, max);
    } catch (IloCplex::Exception& e) {
        cerr << e << endl;
    }
    // SecondaryStructure bestSSO1 = myMOIP.solve_objective(1, -max, max);
    // SecondaryStructure bestSSO2 = myMOIP.solve_objective(2, -max, max);
    // double             bestObj2 = bestSSO2.get_objective_score(2);

    // extend to the whole pareto set
    // myMOIP.add_solution(bestSSO1);
    // myMOIP.extend_pareto(bestObj2, max);

    // print the pareto set
    // cout << "Structure \t Free energy score \t Expected accuracy score" << endl;
    // for (uint i = 0; i < myMOIP.get_n_solutions(); i++) {
    //     cout << myMOIP.solution(i).to_string() << endl;
    // }
    // cout << endl;

    return EXIT_SUCCESS;
}
