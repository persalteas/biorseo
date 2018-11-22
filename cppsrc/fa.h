/*
 * $Id$
 *
 * Copyright (C) 2008-2010 Kengo Sato
 *
 * This file comes from Ipknot.
*/

#include <list>
#include <string>

using std::list;
using std::string;

class Fasta
{
  public:
    Fasta() : name_(), seq_(), str_() {}
    Fasta(const string& name, const string& seq, const string& str = "") : name_(name), seq_(seq), str_(str) {}
    Fasta(const Fasta& fa) : name_(fa.name_), seq_(fa.seq_), str_(fa.str_) {}
    Fasta& operator=(const Fasta& fa)
    {
        if (this != &fa) {
            name_ = fa.name_;
            seq_  = fa.seq_;
            str_  = fa.str_;
        }
        return *this;
    }
    const string& name() const { return name_; }
    const string& seq() const { return seq_; }
    string&       seq() { return seq_; }
    const string& str() const { return str_; }
    unsigned int  size() const { return seq_.size(); }

    static unsigned int load(list<Fasta>& data, const char* file);

  private:
    string name_;
    string seq_;
    string str_;
};
