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
		int nframe; //ICNTRL[1], NFILE
		int tstart; //ICNTRL[2], NPRIV
		int first;
		int delta; //ICNTRL[3], NSAVC
		double deltat;
		int dof;
		float tstep; //ICNTRL[10], DELTA
		bool crystal;
		double pbx;
		double pby;
		double pbz;
		bool fixed;
		unsigned int version;
		int natom;
		std::string endian;
		std::string title;
		std::vector<float> x;
		std::vector<float> y;
		std::vector<float> z;
		std::vector<int> fixinx;

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
		void readFrame(std::ifstream &trjin, unsigned int frame);
		std::string getHeader(){return hdr;};
		int getNFrame();
    void setMolecule(Molecule *molin);
};


#endif
