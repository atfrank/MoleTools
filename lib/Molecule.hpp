//Sean M. Law
#ifndef MOLECULE_H
#define MOLECULE_H

#include "Chain.hpp"

#include <vector>

class Molecule {
  private:
    std::vector<Chain> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;

  public:
    static Molecule* readPDB (std::string ifile, int model=0);
    int writePDB ();
    void addAtom(Atom* atmEntry);
    Atom* getAtom(int element); //This should not be pointer!
    void addChain(Chain chnEntry);
    void addResidue(Residue* resEntry);
    Chain* getChain(int element);
    unsigned int getChnVecSize();
    Chain* getLastChainRef();
    unsigned int getResVecSize();
    Residue* getLastResidueRef();
    unsigned int getAtmVecSize();
    Residue* getResidue(int element);
    Atom* getLastAtomRef();
};

#endif
