#include <iostream>
extern "C"
{
	#include <ViennaRNA/part_func_window.h>
}

#include "rna.h"


using namespace Eigen;
using std::cerr;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

RNA::RNA(void) {}

RNA::RNA(string name, string seq, bool verbose)
: verbose_{verbose}, name_(name), seq_(seq), n_(seq.size()),
  pij_(MatrixXf::Zero(n_, n_))
{
	vector<char> unknown_chars;
	bool         contains_T = false;
	for (char c : seq) {
		if (c == 'T' or c == 't') {
			c          = 'U';
			contains_T = true;
		}
		if (base_type(c) == BASE_N) {
			unknown_chars.push_back(c);
		}
	}
	if (contains_T)
		if (verbose_) cout << "\tWARNING: Thymines automatically replaced by uraciles.";
	if (unknown_chars.size() > 0) {
		if (verbose_) cout << "\tWARNING: Unknown chars in input sequence ignored : ";
		for (char c : unknown_chars)
			if (verbose_) cout << c << " ";
		if (verbose_) cout << endl;
	}
	if (verbose_) cout << "\t> Sequence formatted" << endl;

	// Compute using ViennaRNA
	if (verbose_) cout << "\t> Computing pairing probabilities (ViennaRNA's pfl_fold)..." << endl;
	
	const char      *cseq = seq.c_str();
	int window_size = 100;
	int max_bp_span = 150;
	float cutoff = 1e-6;
	vrna_ep_t* results = vrna_pfl_fold(cseq, window_size, max_bp_span, cutoff);

	if (results != NULL)
	{
		while (results->i != 0  &&  results->j != 0)
		{
			pij_(results->i,results->j) = results->p;
			results++;
		}
	}

	else cout << "NULL result returned by vrna_pfl_fold" << endl;
}





void RNA::print_basepair_p_matrix(float theta) const
{
	cout << endl;
	cout << "\t=== -log10(p(i,j)) for each pair (i,j) of nucleotides: ===" << endl << endl;
	cout << "\t" << seq_ << endl;
	uint i = 0;
	for (uint u = 0; u < pij_.rows(); u++) {
		cout << "\t";
		for (uint v = 0; v < pij_.cols(); v++) {
			if (pij_(u, v) < 5e-10)
				cout << " ";
			else if (pij_(u, v) > theta)
				cout << "\033[0;32m" << int(-log10(pij_(u, v))) << "\033[0m";
			else
				cout << int(-log10(pij_(u, v)));
		}
		cout << seq_[i] << endl;
		i++;
	}
	cout << endl << "\t\033[0;32mgreen\033[0m basepairs are kept as decision variables." << endl << endl;
}

base_t RNA::base_type(char x) const
{
	if (x == 'a' or x == 'A') return BASE_A;
	if (x == 'c' or x == 'C') return BASE_C;
	if (x == 'g' or x == 'G') return BASE_G;
	if (x == 'u' or x == 'U') return BASE_U;
	return BASE_N;
}
