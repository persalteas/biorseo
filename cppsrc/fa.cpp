/*
 * $Id$
 *
 * Copyright (C) 2008-2010 Kengo Sato
 *
 * This file comes from Ipknot.
 */


#include "fa.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <cstring>
#include <cassert>

typedef unsigned int uint;

//static
unsigned int Fasta::load(std::list<Fasta>& data, const char* file){

  std::string line, name, seq, str;
  std::ifstream ifs;
  ifs.open(file, std::ios::in);
  while (std::getline(ifs, line)) {
    if (line[0]=='>') {         // header
      if (!name.empty()) {
        assert(str.size()==0 || seq.size()==str.size());
        data.push_back(Fasta(name, seq, str));
      }

      name=line.substr(1);
      seq.clear();
      str.clear();
      continue;
    } 

    if (std::strchr("()[].?xle ", line[0])==NULL) { // seq
      uint i;
      for (i=0; i!=line.size(); ++i)
        if (!isalpha(line[i])) break;
      seq+=line.substr(0, i);
    } else {
      uint i;
      for (i=0; i!=line.size(); ++i)
        if (std::strchr("()[].?xle ", line[i])==NULL) break;
      str+=line.substr(0, i);
    }
  }
  
  if (!name.empty()) 
    data.push_back(Fasta(name, seq, str));

  return data.size();
}

