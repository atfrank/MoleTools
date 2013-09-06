//Sean M. Law

#ifndef MOLECULE_H
#define MOLECULE_H

#include "Chain.hpp"
#include "LinAlg.hpp"
#include "Eigen/Dense"

#include <vector>
#include <map>

//Base class
class Molecule {
  private:
    std::vector<Chain*> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;
		bool copyFlag; //This molecule is a copy if true
		std::map< std::string, std::vector<bool> > storedSel;
		std::string remarks;
		bool iCodeFlag;

  public:
		Molecule(); //Constructor
		~Molecule();
    static Molecule* readPDB (const std::string ifile, const int model=0, const std::string format="");
		static Molecule* readPDB (const std::string ifile, const std::string format);
//    std::string writePDB (bool selFlag=true, bool print=true, bool chnFlag=false);
		std::string writePDB (bool selFlag, bool print, bool chnFlag);
		std::string writePDB (bool selFlag, bool print);
		std::string writePDB (bool chnFlag);
		std::string writePDB ();
		static Molecule* readMol2 (const std::string ifile, const std::string format="");
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
    unsigned int getNAtomSelected(); //Determining this on the fly is a good safety measure
		void setCopyFlag(bool copyFlagIn=false);
		bool getCopyFlag();
		void storeSel(std::string key="tmp");
		void recallSel(std::string key="tmp");
		void eraseSel(std::string key="tmp");
		void zeroCoor();
		void addRemark(const std::string& remin);
		void clearRemark();
		std::string getRemark();
		bool checkICode();
		void setICodeFlag(bool iCodeFlagIn=false);
		bool getICodeFlag();

    double lsqfit (Molecule *refmol, bool transform=true);
		double rmsd (Molecule *refmol);
		void recenter (Molecule *recmol);
    void translate (const double &dx, const double &dy, const double &dz);
    void translate (const Vector &u);
    void rotate (const double &r1c1, const double &r1c2, const double &r1c3, const double &r2c1, const double &r2c2, const double &r2c3, const double &r3c1, const double &r3c2, const double &r3c3);
    void center (bool selFlag=true);
		void modPseudoCenter();
		void pcasso (std::string dsspin="");

		//Virtual functions
		virtual void format();
};

class MoleculeCHARMM: public Molecule{
	public:
		void format();
};

#endif
