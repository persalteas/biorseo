
#ifndef DEF_RNA
#define DEF_RNA

#include <map>
#include <string>
#include <vector>

using std::pair, std::string, std::vector, std::map;

typedef struct Comp_ {
    pair<uint, uint> pos;
    int              score;
    size_t           k;
    Comp_(pair<int, int> p, int s) : pos(p), score(s) { k = 1 + pos.second - pos.first; }
} Component;
typedef struct {
    string            atlas_id;
    vector<Component> comp;
    bool              reversed;
} Motif;

class RNA
{
  public:
    RNA();
    RNA(string name, string seq);

    int                    get_n_();
    string                 get_name_();
    string                 get_seq_();
    vector<vector<int>>    get_type_();
    int                    get_type(int i, int j);
    int                    get_type(int i);
    vector<pair<int, int>> get_coord_();
    pair<int, int>         get_coord(int i);
    int                    get_coordF(int i);
    int                    get_coordS(int i);
    int                    find_coord(pair<int, int>);
    vector<vector<float>>  get_pij_();
    float                  get_pij(int i, int j);
    float                  get_pij(int i);
    int                    get_err_();
    uint                   get_RNA_length() const;

    bool check_seq(string seq);
    void format();

  private:
    string                          name_;  /*name of the rna*/
    string                          seq_;   /*sequence of the rna*/
    int                             n_;     /*length of the rna*/
    vector<vector<int>>             type_;  /*vector of base pair types*/
    vector<pair<int, int>>          coord_; /*vector of base pair coordinates*/
    vector<vector<float>>           pij_;   /*vector of probabilities*/
    uint                            nBP_;   /*number of possible base pair*/
};

inline int                    RNA::get_n_() { return n_; }
inline string                 RNA::get_name_() { return name_; }
inline string                 RNA::get_seq_() { return seq_; }
inline vector<vector<int>>    RNA::get_type_() { return type_; }
inline int                    RNA::get_type(int i, int j) { return type_[i][j]; }
inline int                    RNA::get_type(int i) { return type_[get_coord(i).first][get_coord(i).second]; }
inline vector<pair<int, int>> RNA::get_coord_() { return coord_; }
inline pair<int, int>         RNA::get_coord(int i) { return coord_[i]; }
inline int                    RNA::get_coordF(int i) { return coord_[i].first; }
inline int                    RNA::get_coordS(int i) { return coord_[i].second; }
inline vector<vector<float>>  RNA::get_pij_() { return pij_; }
inline float                  RNA::get_pij(int i, int j) { return pij_[i][j]; }
inline float                  RNA::get_pij(int i) { return pij_[get_coord(i).first][get_coord(i).second]; }
inline uint                   RNA::get_RNA_length() const { return nBP_; }

#endif
