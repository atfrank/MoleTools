//Sean M. Law
#ifndef MOLECULE_H
#define MOLECULE_H

#include "Chain.hpp"

#include <vector>
#include <Eigen/Dense>

class Molecule {
  private:
    std::vector<Chain*> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;
		bool copyFlag; //This molecule is a copy if true

  public:
		~Molecule();
    static Molecule* readPDB (std::string ifile, int model=0);
    void writePDB (bool selFlag=true);
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
		void setCopyFlag(bool copyFlagIn=false);
		bool getCopyFlag();

    double lsqfit (Molecule *refmol, bool transform=true);
		double rmsd (Molecule *refmol);
		void recenter (Molecule *recmol);
};

#endif
