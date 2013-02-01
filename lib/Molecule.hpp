//Sean M. Law
#ifndef MOLECULE_H
#define MOLECULE_H

#include "Chain.hpp"

#include <vector>

class Molecule {
  private:
    std::vector<Chain> chnVec;
    std::vector<Residue> resVec;
    std::vector<Atom> atmVec;

  public:
  
    int readPDB (std::string ifile, int model=0);
    int writePDB ();
    Atom getAtom(int atmnum);
    unsigned int getAtmVecSize();
};

#endif
