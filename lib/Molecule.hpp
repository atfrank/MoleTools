//Sean M. Law
#ifndef MOLECULE_H
#define MOLECULE_H

#include "Chain.hpp"

#include <vector>

class Molecule {
  private:
    std::vector<Chain*> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;

  public:
		~Molecule();
    static Molecule* readPDB (std::string ifile, int model=0);
    int writePDB ();
    Molecule* clone(bool selFlag=true, bool keep=true);
		Molecule* copy(bool selFlag=true);
    void addAtom(Atom* atmEntry);
    Atom* getAtom(int element); 
    void addChain(Chain* chnEntry);
    void addResidue(Residue* resEntry);
    Chain* getChain(int element);
    unsigned int getChnVecSize();
    unsigned int getResVecSize();
    unsigned int getAtmVecSize();
    std::vector<Atom*> getAtmVec();
    Residue* getResidue(int element);
    void selAll();
    void deselAll();
    void select(std::string sel);
    unsigned int getNAtom();
    unsigned int getNAtomSelected();

    //Analysis Functions
    void lsqfit (Molecule *refmol);
};

#endif
