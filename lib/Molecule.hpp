//Sean M. Law

#ifndef MOLECULE_H
#define MOLECULE_H

#include "Chain.hpp"
#include "LinAlg.hpp"
#include "Eigen/Dense"

#include <vector>
#include <map>

class Molecule {
  private:
    std::vector<Chain*> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;
		bool copyFlag; //This molecule is a copy if true
		std::map< std::string, std::vector<bool> > storedSel;

  public:
		~Molecule();
    static Molecule* readPDB (std::string ifile, int model=0);
    std::string writePDB (bool selFlag=true, bool print=true);
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
		void storeSel(std::string key="tmp");
		void recallSel(std::string key="tmp");
		void eraseSel(std::string key="tmp");
		void zeroCoor();

    double lsqfit (Molecule *refmol, bool transform=true);
		double rmsd (Molecule *refmol);
		void recenter (Molecule *recmol);
    void translate (const double &dx, const double &dy, const double &dz);
    void translate (const Vector &u);
    void rotate (const double &r1c1, const double &r1c2, const double &r1c3, const double &r2c1, const double &r2c2, const double &r2c3, const double &r3c1, const double &r3c2, const double &r3c3);
    void center (bool selFlag=true);
};

#endif
