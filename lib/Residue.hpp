//Sean M. Law

#ifndef RESIDUE_H
#define RESIDUE_H

#include "Atom.hpp"

#include <vector>

class Residue {
  private:
    std::string resname;
    int resnum;
    std::string chainid;
    int start;
    int end;
    std::string segid;
    std::vector<Atom> *atmVec;
};

#endif
