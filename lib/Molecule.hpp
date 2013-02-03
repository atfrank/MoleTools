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
    void addAtom(Atom atmEntry);
    Atom getAtom(int element);
    unsigned int getAtmVecSize();
    void addChain(Chain chnEntry);
    Chain getChain(int element);
    unsigned int getChnVecSize();
    Chain& getLastChainRef();
    unsigned int getResVecSize();
    Residue& getLastResidueRef();
    Atom& getLastAtomRef();
};

#endif
