//Sean M. Law

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "Constants.hpp"
#include "Molecule.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <climits>
#include <cstdlib>

class Trajectory {
  private:
    std::string format;
    bool swab; //Swap bytes
    Molecule *mol;

		std::string hdr;
		int nframe; //ICNTRL[1], Number of frames
		int npriv; //ICNTRL[2], Number of previous integration steps
		int nsavc; //ICNTRL[3], Frequency for saving of frames
		int nstep; //ICNTRL[4],Number of steps in the run that created this file
		int qvelocity; //ICNTRL[5], Velocity flag
		
		int dof; //ICNTRL[8], Degrees of freedom
		int nfixed; //ICNTRL[9], Number of fixed atoms
    float tstep; //ICNTRL[10], AKMA units
		int qcrystal; //ICNTRL[11], Crystal lattice/Periodic boundaries flag
		int q4d; //ICNTRL[12], 4D trajectory flag
		int qcharge; //ICNTRL[13], Fluctuating charges flag
		int qcheck; //ICNTRL[14], Consistency check flag, See "NOCHeck" keyword

		int version; //ICNTRL[20], CHARMM version
	
		std::string title1;
    std::string title2;
		int natom;
		std::vector<int> fixinx;

		double pbx;
		double pby;
		double pbz;
		double pbalpha;
		double pbbeta;
		double pbgamma;
		
		std::vector<float> x;
		std::vector<float> y;
		std::vector<float> z;

		std::string endian;

		typedef union {
			unsigned int ui;
			int i;
			char c[4];
			float f;
		} binbuf;

  public:
    Trajectory ();
    bool findFormat(std::ifstream &trjin);
		template <class BinBuf> 
			BinBuf* readFortran(std::ifstream &trjin, BinBuf *buffer, int &length);
    template <class BinBuf>
      void writeFortran(std::ofstream &trjout, BinBuf *buffer, int &length);
    void writeFortran(std::ofstream &trjout);
		void clearHeader();
		void readHeader(std::ifstream &trjin);
    void writeHeader(std::ofstream &trjout);
		void showHeader();
		void cloneHeader(Trajectory *ftrjin);
		void readFrame(std::ifstream &trjin, unsigned int frame);
		std::string getHeader(){return hdr;};
    void setMolecule(Molecule *molin);

		//Get
    std::string getFormat();
		std::string getHdr();
    int getNFrame();
    int getNPriv(); 
    int getNSavc();
    int getNStep();
    int getQVelocity();

    int getDOF(); 
    int getNFixed(); 
    float getTStepAKMA();
    double getTStepPS();
    int getQCrystal(); 
    int getQ4D(); 
    int getQCharge(); 
    int getQCheck(); 

    int getVersion(); 

    std::string getTitle1();
    std::string getTitle2();
    int getNAtom();
		unsigned int getFixInxVecSize();
		int getFixInx(int element);	

		//Set
    void setHdr(const std::string &hdrin);
    void setNFrame(const int &nframein);
    void setNPriv(const int &nprivin);
    void setNSavc(const int &nsavcin);
    void setNStep(const int &nstepin);
    void setQVelocity(const int &qvelocityin);

    void setDOF(const int &dofin);
    void setNFixed(const int &nfixedin);
    void setTStep(const double &tstepin); //Picosecond input units
    void setTStep(const float &tstepin); //AKMA input units
    void setQCrystal(const int &qcrystalin);
    void setQ4D(const int &q4din);
    void setQCharge(const int &qchargein);
    void setQCheck(const int &qcheckin);

    void setVersion(const int &versionin);

    void setTitle1(const std::string &title1in);
    void setTitle2(const std::string &title2in);
    void setNAtom(const int &natomin);
		void addFixInx(const int &elementin);
		void clearFixInx();
};


#endif
