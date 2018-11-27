
#ifndef DEF_RNA
#define DEF_RNA

#include <map>
#include <sstream>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::string;
using std::vector;

typedef struct Comp_ {
    pair<uint, uint> pos;
    int    score;
    size_t k;
    Comp_(pair<int, int> p, int s) : pos(p), score(s) { k = 1 + pos.second - pos.first; }
} Component;
typedef struct {
    string            atlas_id;
    vector<Component> comp;
    bool              reversed;
    string            pos_string(void) const
    {
        std::stringstream s;
        for (auto c : comp) {
            s << c.pos.first << '-' << c.pos.second << ' ';
        }
        return s.str();
    }
} Motif;

class RNA
{
    public:
    RNA();
    RNA(string name, string seq);

    int                   get_n();
    string                get_name();
    string                get_seq();
    vector<vector<float>> get_pij();
    float get_pij(int i, int j);
    float get_pij(int i);
    int  get_err();
    uint get_RNA_length() const;
    void print_basepair_p_matrix(float theta) const;

    bool check_seq(string seq);
    void format();

    private:
    string                name_; /* name of the rna */
    string                seq_;  /* sequence of the rna */
    int                   n_;    /* length of the rna */
    vector<vector<float>> pij_;  /* vector of probabilities */
};

inline int                   RNA::get_n() { return n_; }
inline string                RNA::get_name() { return name_; }
inline string                RNA::get_seq() { return seq_; }
inline vector<vector<float>> RNA::get_pij() { return pij_; }
inline float RNA::get_pij(int i, int j) { return pij_[i][j]; }
inline uint RNA::get_RNA_length() const { return n_; }

#endif
