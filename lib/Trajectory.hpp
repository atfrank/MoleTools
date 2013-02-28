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
		bool qvelocity; //ICNTRL[5], Velocity flag
		
		int dof; //ICNTRL[8], Degrees of freedom
		int nfixed; //ICNTRL[9], Number of fixed atoms
		double tstep; //ICNTRL[10], AKMA units 
		bool qcrystal; //ICNTRL[11], Crystal lattice/Periodic boundaries flag
		bool q4d; //ICNTRL[12], 4D trajectory flag
		bool qcharge; //ICNTRL[13], Fluctuating charges flag
		bool qcheck; //ICNTRL[14], Consistency check flag, See "NOCHeck" keyword

		int version; //ICNTRL[20], CHARMM version
	
		std::string title;
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
		void clearHeader();
		void readHeader(std::ifstream &trjin);
		void showHeader();
		void cloneHeader(Trajectory *ftrjin);
		void readFrame(std::ifstream &trjin, unsigned int frame);
		std::string getHeader(){return hdr;};
    void setMolecule(Molecule *molin);

		//Get
		std::string getHdr();
    int getNFrame();
    int getNPriv(); 
    int getNSavc();
    int getNStep();
    bool getQVelocity();

    int getDOF(); 
    int getNFixed(); 
    double getTStep();  
    bool getQCrystal(); 
    bool getQ4D(); 
    bool getQCharge(); 
    bool getQCheck(); 

    int getVersion(); 

    std::string getTitle();
    int getNAtom();
		unsigned int getFixInxVecSize();
		int getFixInx(int element);		
		//Set
};


#endif
